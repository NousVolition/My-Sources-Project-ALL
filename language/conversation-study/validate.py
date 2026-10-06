"""Check provenance and evidence links, not the truth of interpretations."""
import hashlib
import json
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_bundle(root=HERE):
    names = ('manifest', 'evidence', 'word-paths', 'confusion-tracks', 'record-index')
    return {name: json.loads((root / f'{name}.json').read_text(encoding='utf-8')) for name in names}


def unique(records, label):
    ids = [row['id'] for row in records]
    require(len(ids) == len(set(ids)), f'Duplicate {label} ID')
    return set(ids)


def validate_data(data, source_bytes):
    manifest = data['manifest']
    sources = {}
    unique(manifest['sources'], 'source')
    for source in manifest['sources']:
        sid = source['id']
        raw = source_bytes[sid]
        require(hashlib.sha256(raw).hexdigest() == source['sha256'], f'Source hash mismatch: {sid}')
        require(len(raw) == source['bytes'], f'Source byte count mismatch: {sid}')
        sources[sid] = raw.decode('utf-8-sig').splitlines()
        require(len(sources[sid]) == source['lines'], f'Source line count mismatch: {sid}')

    def span(ref):
        sid, start, end = ref
        require(sid in sources, f'Unknown source: {sid}')
        require(isinstance(start, int) and isinstance(end, int) and 1 <= start <= end <= len(sources[sid]), f'Invalid source span: {ref}')
        return '\n'.join(sources[sid][start - 1:end])

    def walk_refs(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key == 'ref':
                    span(item)
                elif key == 'refs' or key.endswith('_refs'):
                    for reference in item:
                        span(reference)
                else:
                    walk_refs(item)
        elif isinstance(value, list):
            for item in value:
                walk_refs(item)

    for overlap in manifest['overlap']:
        old = [overlap['previous_source'], *overlap['previous_lines']]
        new = [overlap['new_source'], *overlap['new_lines']]
        require(span(old) == span(new), 'Recorded overlap is not exact')
    evidence = data['evidence']
    categories = unique(evidence['categories'], 'category')
    episodes = unique(evidence['events'], 'episode')
    for event in evidence['events']:
        require(set(event['tags']) <= categories, f'Unknown category: {event["id"]}')
        require(bool(event['refs']), f'Missing evidence: {event["id"]}')
        walk_refs(event)
    paths = data['word-paths']['paths']
    unique(paths, 'path')
    fragments = 0
    for path in paths:
        for node in path['nodes']:
            require(bool(node['quote']) and node['quote'] in span(node['ref']), f'Quote mismatch: {path["id"]}')
            fragments += 1
    context = data['confusion-tracks']
    unique(context['tracks'], 'context track')
    walk_refs(context)
    for track in context['tracks']:
        require(set(track['related_episodes']) <= episodes, f'Unknown episode: {track["id"]}')
        for step in track['steps']:
            for quote in step['quotes']:
                require(bool(quote['text']) and quote['text'] in span(quote['ref']), f'Quote mismatch: {track["id"]}')
                fragments += 1
    for reference in context['same_request']['refs']:
        require(span(reference) == context['same_request']['text'], 'Repeated request mismatch')
    # Both answer premises must also match; the request alone would be a weaker comparison.
    for first, second in ((1901, 2213), (1902, 2214)):
        a = span(['source-03', first, first]).split('"')
        b = span(['source-03', second, second]).split('"')
        require(len(a) >= 3 and len(b) >= 3 and a[1] == b[1], 'Repeated comparison premise mismatch')
    index = data['record-index']
    records = unique(index['records'], 'indexed record')
    for record in index['records']:
        span([record['source'], record['user_start'], record['user_end']])
        span([record['source'], record['reply_start'], record['reply_end']])
    for alignment in index['reconstruction_alignment']:
        require(alignment['candidate_record'] is None or alignment['candidate_record'] in records, 'Unknown reconstruction candidate')
        span(['source-02', alignment['line'], alignment['line']])
    summary = manifest['index_summary']
    require(summary['selected_evidence_episodes'] == len(episodes), 'Episode total mismatch')
    require(summary['selected_word_paths'] == len(paths), 'Path total mismatch')
    require(summary['context_comparisons_reusing_existing_evidence'] == len(context['tracks']), 'Context total mismatch')
    require(summary['earlier_exchange_blocks'] == len(records), 'Index total mismatch')
    require(summary['reconstruction_entries_including_one_empty_truncation'] == len(index['reconstruction_alignment']), 'Reconstruction total mismatch')
    require(summary['question_alignment'] == dict(Counter(a['question_match'] for a in index['reconstruction_alignment'])), 'Alignment totals mismatch')
    return {'sources': len(sources), 'observations': len(episodes), 'categories': len(categories), 'paths': len(paths), 'context_tracks': len(context['tracks']), 'exact_fragments': fragments}


class Page(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.links = []
        self.ids = set()
        self.text = []
        self.source_lines = []
        self.source_span = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            require(attrs['id'] not in self.ids, f'Duplicate HTML ID: {attrs["id"]}')
            self.ids.add(attrs['id'])
        if tag == 'a' and 'href' in attrs:
            self.links.append(attrs['href'])
        if tag == 'span' and attrs.get('class') == 'text':
            self.source_span = []

    def handle_data(self, value):
        self.text.append(value)
        if self.source_span is not None:
            self.source_span.append(value)

    def handle_endtag(self, tag):
        if tag == 'span' and self.source_span is not None:
            self.source_lines.append(''.join(self.source_span))
            self.source_span = None


def validate_pages(root=HERE):
    root = root.resolve()
    pages = {}
    for path in root.rglob('*.html'):
        page = Page()
        page.feed(path.read_text(encoding='utf-8'))
        pages[path.resolve()] = page
    count = 0
    for path, page in pages.items():
        for href in page.links:
            url = urlsplit(href)
            if url.scheme in ('http', 'https'):
                continue
            target = (path.parent / unquote(url.path)).resolve() if url.path else path
            require(target.is_relative_to(root) and target.is_file(), f'Broken local link: {href}')
            if url.fragment:
                require(target in pages and unquote(url.fragment) in pages[target].ids, f'Broken anchor: {href}')
            count += 1
    data = read_bundle(root)
    for source in data['manifest']['sources']:
        sid = source['id']
        expected = (root / 'sources' / f'{sid}.txt').read_text(encoding='utf-8-sig').splitlines()
        require(pages[root / 'sources' / f'{sid}.html'].source_lines == expected, f'Source display differs: {sid}')
    page = pages[root / 'index.html']
    plain = ''.join(page.text)
    for event in data['evidence']['events']:
        require(event['id'] in page.ids and event['observation'] in plain and event['limit'] in plain, f'Interactive observation differs: {event["id"]}')
    return count


def main():
    data = read_bundle()
    raw = {s['id']: (HERE / 'sources' / f'{s["id"]}.txt').read_bytes() for s in data['manifest']['sources']}
    result = validate_data(data, raw)
    result['local_html_links'] = validate_pages()
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
