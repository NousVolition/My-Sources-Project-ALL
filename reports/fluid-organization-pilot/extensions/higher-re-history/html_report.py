import html,re
def inline(s):
 s=html.escape(s);s=re.sub(r'`([^`]+)`',r'<code>\1</code>',s);s=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',s);return re.sub(r'\[([^]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',s)
def render(text):
 out=[];lines=text.splitlines();i=0
 while i<len(lines):
  line=lines[i]
  if not line.strip():i+=1;continue
  if line.startswith('```'):
   code=[];i+=1
   while i<len(lines) and not lines[i].startswith('```'):code.append(lines[i]);i+=1
   out.append('<pre><code>'+html.escape('\n'.join(code))+'</code></pre>');i+=1;continue
  if line.startswith('|'):
   rows=[]
   while i<len(lines) and lines[i].startswith('|'):rows.append([c.strip() for c in lines[i].strip('|').split('|')]);i+=1
   out.append('<div class="table"><table><thead><tr>'+''.join('<th>'+inline(c)+'</th>' for c in rows[0])+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+inline(c)+'</td>' for c in row)+'</tr>' for row in rows[2:])+'</tbody></table></div>');continue
  if line.startswith('!['):
   m=re.match(r'!\[([^]]*)\]\(([^)]+)\)',line);out.append('<figure><img alt="'+html.escape(m[1])+'" src="'+html.escape(m[2])+'"></figure>');i+=1;continue
  if line.startswith('#'):
   level=len(line)-len(line.lstrip('#'));out.append(f'<h{level}>'+inline(line[level:].strip())+f'</h{level}>');i+=1;continue
  if line.startswith('- '):
   items=[]
   while i<len(lines) and lines[i].startswith('- '):items.append('<li>'+inline(lines[i][2:])+'</li>');i+=1
   out.append('<ul>'+''.join(items)+'</ul>');continue
  para=[line];i+=1
  while i<len(lines) and lines[i].strip() and not lines[i].startswith(('#','|','![','- ','```')):para.append(lines[i]);i+=1
  out.append('<p>'+inline(' '.join(para))+'</p>')
 return '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>History prediction at increased Reynolds number</title><style>body{font:16px/1.65 system-ui,Segoe UI,sans-serif;color:#24323b;max-width:1120px;margin:auto;padding:32px}h1{line-height:1.2;font-size:2.15rem}h2{margin-top:2.1em;line-height:1.3}a{color:#086580}code{background:#eef3f6;padding:2px 4px}pre{background:#eef3f6;padding:16px;overflow:auto;font-size:13px}.table{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.91em}th,td{border-bottom:1px solid #d7e0e5;text-align:left;padding:8px 10px}th{background:#edf3f6}img{width:100%;height:auto}figure{margin:25px 0}</style><main>'+''.join(out)+'</main></html>'
