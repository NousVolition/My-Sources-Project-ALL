"""Render GitHub reading pages from the saved annotations. Standard library only."""
import argparse
import html
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE / name).read_text(encoding='utf-8'))


def text(value):
    return html.escape(str(value), quote=False).replace('|', '&#124;')


def ref(span):
    source, start, end = span
    label = f'{source} L{start}' + (f'–{end}' if end != start else '')
    return f'[{label}](sources/{source}.txt#L{start}' + (f'-L{end}' if end != start else '') + ')'


def refs(spans):
    return ' · '.join(ref(span) for span in spans)


def blockquote(value):
    return '\n'.join('> ' + text(line) for line in value.splitlines())


def heading(title, intro):
    return [f'# {title}', '[← Study home](README.md)', intro]


def render():
    evidence = load('evidence.json')
    context = load('confusion-tracks.json')
    paths = load('word-paths.json')
    manifest = load('manifest.json')
    categories = {c['id']: c['name'] for c in evidence['categories']}
    pages = {}
    lines = heading('Confusion and context', context['method'] + '\n\n' + context['causal_limit'])
    lines.append(' · '.join(f'[{t["id"]}](#{t["id"].lower()})' for t in context['tracks']))
    for track in context['tracks']:
        lines.extend([f'## {track["id"]}', f'### {text(track["title"])}', text(track['summary'])])
        if track['id'] == 'C01':
            lines.extend(['**The identical request, preserved verbatim:**', blockquote(context['same_request']['text']), refs(context['same_request']['refs'])])
        for step in track['steps']:
            lines.append(f'#### {text(step["label"])}')
            for label, key, rk in [('Before', 'before', 'before_refs'), ('Request', 'request', 'request_refs'), ('Answer', 'answer', 'answer_refs'), ('Next clarification', 'next_clarification', 'next_refs')]:
                lines.append(f'**{label}:** {text(step[key])}\n\n{refs(step[rk])}')
            lines.append(f'**Visible status:** {text(step["status"])}')
            for quote in step['quotes']:
                lines.extend([blockquote(quote['text']), ref(quote['ref'])])
        for label, key in [('Changed', 'what_changed'), ('Carried forward', 'what_stayed'), ('Outcome', 'outcome'), ('Reading limit', 'limit')]:
            lines.append(f'**{label}:** {text(track[key])}')
        lines.append('Related observations: ' + ' · '.join(f'[{eid}](OBSERVATIONS.md#{eid.lower()})' for eid in track['related_episodes']))
    pages['CONTEXT.md'] = '\n\n'.join(lines) + '\n'
    lines = heading('Word paths', paths['scope'])
    for path in paths['paths']:
        lines.extend([f'## {path["id"]}', f'### {text(path["title"])}', text(path['reading'])])
        for node in path['nodes']:
            lines.extend([f'**{text(node["speaker"])} · {text(node["label"])}**', blockquote(node['quote']), ref(node['ref'])])
        lines.append('**Reading limit:** ' + text(path['limit']))
    pages['PATHS.md'] = '\n\n'.join(lines) + '\n'
    lines = heading('Source-linked observations', evidence['scope'])
    for event in evidence['events']:
        lines.extend([f'## {event["id"]}', f'### {text(event["title"])}', '**Categories:** ' + ' · '.join(categories[tag] for tag in event['tags']), '**Observed:** ' + text(event['observation']), '**Reading limit:** ' + text(event['limit']), refs(event['refs'])])
        if event.get('external'):
            lines.append(f'[Citation identifier supplied for inspection]({event["external"]})')
    pages['OBSERVATIONS.md'] = '\n\n'.join(lines) + '\n'
    lines = heading('Source snapshots', manifest['completeness'])
    lines.append(manifest['publication'])
    for source in manifest['sources']:
        sid = source['id']
        lines.extend([f'## {sid}', f'[{source["lines"]:,} lines of unchanged supplied text](sources/{sid}.txt) · {source["bytes"]:,} bytes', source['role'], f'SHA-256: `{source["sha256"]}`'])
    lines.extend(['## Known overlap', 'Source-03 lines 1–565 exactly repeat source-02 lines 777–1341. Keep both originals; do not count the repeated span as new independent exchanges. The continuation starts at [source-03 line 566](sources/source-03.txt#L566).', '## Earlier reconstruction alignment', 'The index contains 88 earlier exchange blocks and 134 embedded reconstruction entries, including one empty/truncated entry. Of the candidate user-text alignments, 128 are exact, three near, two partial, and one unmatched. A reconstruction remains an AI-authored reconstruction even when text matches.', '[Manifest and complete rules](manifest.json) · [Candidate reconstruction alignments](record-index.json)', '## Reading rules'])
    lines.extend(f'- **{text(key.title())}:** {text(value)}' for key, value in manifest['rules'].items())
    pages['SOURCES.md'] = '\n\n'.join(lines) + '\n'
    lines = heading('Categories without the reply text', 'Columns are selected episodes, not individual replies. A filled cell records a provisional tag. A blank means no tag was assigned. Four later categories have not been applied retrospectively; blanks are not evidence of absence.')
    lines.append('![Category map of 53 selected episodes across 13 overlapping categories](category-map.svg)')
    lines.extend(['## Definitions', '| Category | Reading definition |', '| --- | --- |'])
    lines[-2:] = ['\n'.join(lines[-2:])]
    lines.append('\n'.join(f'| {text(c["name"])} | {text(c["definition"])} |' for c in evidence['categories']))
    # Keep the complete table contiguous for GitHub Markdown.
    lines[-2:] = ['\n'.join(lines[-2:])]
    lines.append('## Inspect a section')
    for title, events in [('Earlier observations · E01–E25', evidence['events'][:25]), ('Later observations · E26–E53', evidence['events'][25:])]:
        lines.append(f'### {title}')
        lines.append(' · '.join(f'[{e["id"]}](OBSERVATIONS.md#{e["id"].lower()})' for e in events))
    lines.append('This graphic is a view of `evidence.json`, not another dataset. Tags and episode order are unchanged. [Read the evidence](OBSERVATIONS.md).')
    pages['CATEGORIES.md'] = '\n\n'.join(lines) + '\n'
    # Four compact panels retain legible labels on GitHub and mobile.
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="980" height="1490" viewBox="0 0 980 1490" role="img" aria-labelledby="title desc">', '<title id="title">Conversation category map</title><desc id="desc">Thirteen rows of categories across fifty-three selected episodes, split into four panels. A dot indicates an assigned tag; an empty cell is unassigned.</desc>', '<rect width="980" height="1490" fill="#f6f8fc"/>', '<g font-family="Arial, sans-serif" fill="#182c42">', '<text x="30" y="40" font-size="24" font-weight="bold">Categories through the conversation</text>', '<text x="30" y="67" font-size="15">Selected episodes • overlapping tags • empty cells are unassigned</text>']
    groups = [evidence['events'][:13], evidence['events'][13:25], evidence['events'][25:39], evidence['events'][39:]]
    for panel, events in enumerate(groups):
        top = 110 + panel * 338
        svg.append(f'<text x="30" y="{top}" font-size="16" font-weight="bold">{events[0]["id"]}–{events[-1]["id"]}</text>')
        for col, event in enumerate(events):
            x = 276 + col * 46
            svg.append(f'<text x="{x}" y="{top + 27}" font-size="12" text-anchor="middle">{event["id"]}</text>')
        for row, cat in enumerate(evidence['categories']):
            y = top + 48 + row * 20
            svg.append(f'<text x="30" y="{y + 4}" font-size="13">{html.escape(cat["name"])}</text>')
            for col, event in enumerate(events):
                x = 276 + col * 46
                svg.append(f'<rect x="{x - 18}" y="{y - 8}" width="36" height="16" rx="3" fill="#e4eaf1"/>')
                if cat['id'] in event['tags']:
                    svg.append(f'<circle cx="{x}" cy="{y}" r="5" fill="#167a78"><title>{event["id"]}: {html.escape(cat["name"])}</title></circle>')
    svg.append('<text x="30" y="1470" font-size="13">Four categories were added later. Earlier episodes have not been rescored for them.</text></g></svg>')
    pages['category-map.svg'] = '\n'.join(svg) + '\n'
    return pages


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Fail if checked-in reading pages are stale.')
    args = parser.parse_args()
    stale = []
    for name, contents in render().items():
        path = HERE / name
        if args.check:
            if not path.exists() or path.read_bytes() != contents.encode('utf-8'):
                stale.append(name)
        else:
            path.write_text(contents, encoding='utf-8', newline='\n')
    if stale:
        raise SystemExit('Rebuild stale views: ' + ', '.join(stale))
    print('GitHub reading pages ' + ('match saved data.' if args.check else 'rebuilt.'))


if __name__ == '__main__':
    main()
