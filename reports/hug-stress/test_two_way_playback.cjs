// Requires Node, Playwright and Chrome. Optional argument: rendered iframe preview.
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const {pathToFileURL}=require('node:url');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
let browser;
(async()=>{
 browser=await chromium.launch({channel:'chrome',headless:true});
 const page=await browser.newPage({viewport:{width:800,height:1000}}),errors=[],layouts=[];
 page.on('pageerror',e=>errors.push(String(e)));
 const target=path.resolve(process.argv[2]||path.join(__dirname,'two-way-report.html'));
 await page.goto(pathToFileURL(target).href);
 let frame=page.mainFrame();const inline=Boolean(await page.locator('iframe').count());
 if(inline){await page.frameLocator('iframe').locator('#hug-two-way').waitFor();frame=page.frames().find(f=>f.parentFrame());}
 const data=await frame.evaluate(()=>JSON.parse(document.getElementById('ht-data').textContent));
 const snap=()=>frame.evaluate(()=>document.getElementById('hug-two-way').hugSnapshot());
 let selections=0;
 for(const c of data.cases){
  await frame.locator('#ht-case').selectOption(c.key);
  const last=data.records[c.one].samples.length-1;
  for(const index of [0,Math.floor(last/2),last]){
   await frame.locator('#ht-time').evaluate((el,i)=>{el.value=i;el.dispatchEvent(new Event('input'));},index);
   const s=await snap();assert.equal(s.index,index);assert.equal(s.key,c.key);
   assert.deepEqual(s.samples,[data.records[c.one].samples[index],data.records[c.two].samples[index]]);
   assert(s.arms.flatMap(a=>[...a.left.flat(),...a.right.flat()]).every(Number.isFinite));selections++;
  }
 }
 await frame.locator('#ht-case').selectOption('closed-0-8');
 for(const width of [320,430,736,1200]){
  await page.setViewportSize({width,height:1000});await page.waitForTimeout(50);
  const layout=await frame.evaluate(()=>({width:innerWidth,scroll:document.documentElement.scrollWidth}));
  assert(layout.scroll<=layout.width,JSON.stringify(layout));layouts.push(layout);
 }
 await page.setViewportSize({width:800,height:1000});await page.waitForTimeout(200);
 const visible=await frame.evaluate(()=>[...document.querySelectorAll('#hug-two-way svg, #ht-play, #ht-case')].map(el=>({name:el.id,left:el.getBoundingClientRect().left,right:el.getBoundingClientRect().right,width:innerWidth})));
 assert(visible.every(x=>x.left>=0&&x.right<=x.width),JSON.stringify(visible));
 await page.screenshot({path:path.join(__dirname,'two-way-data/playback-preview.png')});
 await frame.locator('#ht-play').click();await page.waitForTimeout(150);
 assert((await snap()).index>0);await frame.locator('#ht-play').click();assert.equal((await snap()).running,false);
 await frame.evaluate(()=>{
  const queue=new Map();let id=0;
  window.requestAnimationFrame=f=>{queue.set(++id,f);return id;};window.cancelAnimationFrame=id=>queue.delete(id);
  window.__step=t=>{const fs=[...queue.values()];queue.clear();for(const f of fs)f(t);};
  Object.defineProperty(performance,'now',{value:()=>1000000});
 });
 await frame.locator('#ht-time').evaluate(el=>{el.value=0;el.dispatchEvent(new Event('input'));});
 await frame.locator('#ht-play').click();await frame.evaluate(()=>{__step(0);__step(50);__step(-1);});
 assert.equal((await snap()).index,2);
 await frame.evaluate(()=>__step(100000));assert.equal((await snap()).index,800);assert.equal((await snap()).running,false);
 for(const index of ['bad',{},[],true,NaN,Infinity,-10,99999,15.6]){
  await frame.evaluate(index=>window.dispatchEvent(new CustomEvent('openai:set_globals',{detail:{globals:{widgetState:{privateContent:{schema:'two-way-hug-v1',key:'closed-0-8',index}}}}})),index);
  const expected=Number.isFinite(index)?Math.max(0,Math.min(800,Math.round(index))):0;
  assert.equal((await snap()).index,expected);
 }
 await page.emulateMedia({reducedMotion:'reduce'});await frame.locator('#ht-play').click();assert.equal((await snap()).index,800);
 await page.emulateMedia({colorScheme:'dark'});await page.screenshot({path:path.join(__dirname,'two-way-data/playback-preview-dark.png')});
 assert.deepEqual(errors,[]);
 const result={surface:inline?'sandboxed iframe':'offline report',comparisons:data.cases.length,sampleSelections:selections,clockMismatch:true,invalidSavedState:true,reducedMotion:true,layouts,errors};
 fs.writeFileSync(path.join(__dirname,inline?'two-way-data/browser-verification-inline.json':'two-way-data/browser-verification.json'),JSON.stringify(result,null,2)+'\n');
 console.log(JSON.stringify(result));await browser.close();
})().catch(async error=>{console.error(error);if(browser)await browser.close();process.exitCode=1;});
