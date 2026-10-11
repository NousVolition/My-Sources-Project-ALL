// Node + Playwright + Chrome. Optional argument: sandboxed preview wrapper.
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const {pathToFileURL}=require('node:url');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
let browser;
(async()=>{
 browser=await chromium.launch({channel:'chrome',headless:true});
 const page=await browser.newPage({viewport:{width:800,height:1050}}),errors=[],layouts=[];
 page.on('pageerror',e=>errors.push(String(e)));
 await page.goto(pathToFileURL(path.resolve(process.argv[2]||path.join(__dirname,'report.html'))).href);
 let frame=page.mainFrame();const inline=Boolean(await page.locator('iframe').count());
 if(inline){await page.frameLocator('iframe').locator('#moldable-hug').waitFor();frame=page.frames().find(f=>f.parentFrame());}
 const data=await frame.evaluate(()=>JSON.parse(document.getElementById('mh-data').textContent));
 const snap=()=>frame.evaluate(()=>document.getElementById('moldable-hug').moldingSnapshot());
 let selections=0;
 for(const c of data.cases){
  await frame.locator('#mh-case').selectOption(c.key);
  const last=data.records[c.one].frames.length-1;
  for(const index of [0,Math.floor(last/2),last]){
   await frame.locator('#mh-time').evaluate((el,i)=>{el.value=i;el.dispatchEvent(new Event('input'));},index);
   const s=await snap();assert.equal(s.index,index);assert.equal(s.key,c.key);
   assert.deepEqual(s.frames,[data.records[c.one].frames[index],data.records[c.two].frames[index]]);
   assert(s.frames.flatMap(f=>f.positions.flat()).every(Number.isFinite));
   assert.equal(await frame.locator('#mh-one path').count(),data.triangles.length+2);selections++;
  }
 }
 await frame.locator('#mh-case').selectOption('ellipse');
 for(const width of [320,430,736,1200]){
  await page.setViewportSize({width,height:1050});await page.waitForTimeout(75);
  const layout=await frame.evaluate(()=>({width:innerWidth,scroll:document.documentElement.scrollWidth}));
  assert(layout.scroll<=layout.width,JSON.stringify(layout));layouts.push(layout);
  const labels=await frame.evaluate(()=>[...document.querySelectorAll('#moldable-hug svg text')].map(el=>{
   const r=el.getBoundingClientRect(),p=el.closest('svg').getBoundingClientRect();return {text:el.textContent,contained:r.left>=p.left-1&&r.right<=p.right+1&&r.top>=p.top-1&&r.bottom<=p.bottom+1};}));
  assert(labels.every(x=>x.contained),JSON.stringify(labels.filter(x=>!x.contained)));
 }
 await page.setViewportSize({width:800,height:1050});await page.waitForTimeout(150);
 await page.screenshot({path:path.join(__dirname,'data/playback-preview.png')});
 await frame.locator('#mh-time').evaluate(el=>{el.value=0;el.dispatchEvent(new Event('input'));});
 await frame.locator('#mh-play').click();await page.waitForTimeout(350);assert((await snap()).index>0);
 await frame.locator('#mh-play').click();assert.equal((await snap()).running,false);
 await frame.evaluate(()=>{
  const queue=new Map();let id=0;
  window.requestAnimationFrame=f=>{queue.set(++id,f);return id;};window.cancelAnimationFrame=id=>queue.delete(id);
  window.__step=t=>{const callbacks=[...queue.values()];queue.clear();for(const f of callbacks)f(t);};
  Object.defineProperty(performance,'now',{value:()=>1000000});
 });
 await frame.locator('#mh-time').evaluate(el=>{el.value=0;el.dispatchEvent(new Event('input'));});
 await frame.locator('#mh-play').click();await frame.evaluate(()=>{__step(0);__step(300);__step(-1);});assert.equal((await snap()).index,2);
 const last=data.records[data.cases.find(c=>c.key==='ellipse').one].frames.length-1;
 await frame.evaluate(()=>__step(100000));assert.equal((await snap()).index,last);assert.equal((await snap()).running,false);
 for(const index of ['bad',{},[],true,NaN,Infinity,-10,99999,15.6]){
  await frame.evaluate(index=>window.dispatchEvent(new CustomEvent('openai:set_globals',{detail:{globals:{widgetState:{privateContent:{schema:'moldable-hug-v1',key:'ellipse',index}}}}})),index);
  const expected=Number.isFinite(index)?Math.max(0,Math.min(last,Math.round(index))):0;
  assert.equal((await snap()).index,expected);
 }
 await page.emulateMedia({reducedMotion:'reduce'});await frame.locator('#mh-play').click();assert.equal((await snap()).index,last);
 await page.emulateMedia({colorScheme:'dark'});await page.screenshot({path:path.join(__dirname,'data/playback-preview-dark.png')});
 assert.deepEqual(errors,[]);
 const result={surface:inline?'sandboxed iframe':'offline report',comparisons:data.cases.length,selections,clockMismatch:true,invalidSavedState:true,reducedMotion:true,layouts,errors};
 fs.writeFileSync(path.join(__dirname,inline?'data/browser-verification-inline.json':'data/browser-verification.json'),JSON.stringify(result,null,2)+'\n');
 console.log(JSON.stringify(result));await browser.close();
})().catch(async error=>{console.error(error);if(browser)await browser.close();process.exitCode=1;});
