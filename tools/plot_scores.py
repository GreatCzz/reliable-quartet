import re, os, datetime
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
for f in ['/home/great_czz/.fonts/simsun.ttc']:
    try: font_manager.fontManager.addfont(f)
    except Exception: pass
matplotlib.rcParams['font.family'] = ['SimSun', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

LOG='/tmp/logs'
TZ=datetime.timezone(datetime.timedelta(hours=8))
NOW=datetime.datetime.now(TZ)

def parse_dt(s):
    s=s.replace(' CST ',' ').strip()
    return datetime.datetime.strptime(s,'%a %b %d %H:%M:%S %Y').replace(tzinfo=TZ)

def get_start(out_path):
    if not os.path.exists(out_path): return None
    txt=open(out_path,errors='ignore').read()
    m=re.search(r'start=(\w{3} \w{3} \d+ \d\d:\d\d:\d\d \w+ \d{4})',txt)
    if m: return parse_dt(m.group(1))
    m=re.search(r'=== [\w ]*?cycle \d+ start (\w{3} \w{3} \d+ \d\d:\d\d:\d\d \w+ \d{4})',txt)
    return parse_dt(m.group(1)) if m else None

def events(err_path):
    ev=[]
    for line in open(err_path,errors='ignore'):
        m=re.search(r'\[solver\] sweep \d+:.*?total=(\d+)',line)
        if m: ev.append((int(m.group(1)),225)); continue
        m=re.search(r'\[spr\] round \d+:.*?total=(\d+)',line)
        if m: ev.append((int(m.group(1)),1000)); continue
        m=re.search(r'\[(?:solver|spr)\] initial score = (\d+)',line)
        if m: ev.append((int(m.group(1)),40))
    return ev

runs=[('smoke-7488','smoke 试跑'),('chain-7489','NNI 链'),
      ('spr-7490','SPR 链'),('alt-7491','链 A（NNI 优先）'),
      ('alt2-7492','链 B（SPR 优先）'),('kick-7504','扰动 v1'),('kick2-7505','扰动 v2')]
series={}
for name,label in runs:
    e=os.path.join(LOG,name+'.err'); o=os.path.join(LOG,name+'.out')
    if not os.path.exists(e): continue
    st=get_start(o)
    if st is None: continue
    ev=events(e)
    if not ev: continue
    tot=sum(d for _,d in ev)
    end=min(NOW, datetime.datetime.fromtimestamp(os.path.getmtime(e),TZ))
    span=(end-st).total_seconds()
    scale=span/tot if tot>0 else 1.0
    t=0.0; pts=[]
    for score,d in ev:
        t+=d*scale
        pts.append((st+datetime.timedelta(seconds=t),score))
    series[label]=pts

merged=sorted([p for pts in series.values() for p in pts],key=lambda x:x[0])
best=0; env=[]
for t,s in merged:
    if s>best: best=s; env.append((t,s))

fig,ax=plt.subplots(figsize=(13,6.5))
colors=plt.cm.tab10.colors
for i,(label,pts) in enumerate(series.items()):
    ax.plot([p[0] for p in pts],[p[1]/1e15 for p in pts],'.',ms=2.5,alpha=0.3,
            color=colors[i%10],label=label)
ax.step([p[0] for p in env],[p[1]/1e15 for p in env],where='post',lw=2.6,color='black',label='历史最优包络')
ax.annotate(f'{env[-1][1]:,}',xy=(env[-1][0],env[-1][1]/1e15),xytext=(-230,-30),
            textcoords='offset points',fontsize=9,
            arrowprops=dict(arrowstyle='->',lw=0.7,color='gray'))
for lvl,lab in [(1.65,'1.65e15'),(1.70,'1.70e15')]:
    ax.axhline(lvl,ls=':',lw=0.8,color='gray')
    ax.text(env[0][0],lvl+0.001,lab,fontsize=7.5,color='gray')
ax.set_ylim(1.598,1.746)
ax.axhline(1.6,ls='--',lw=0.9,color='gray'); ax.text(env[0][0],1.601,'1.6e15',fontsize=8,color='gray')
ax.set_xlabel('真实时间（北京时间）'); ax.set_ylabel('分数（×10¹⁵）')
ax.set_title('rep01 分数上升曲线（1000 棵基因树口径，集群任务）')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M',tz=TZ))
ax.xaxis.set_major_locator(mdates.HourLocator(interval=3,tz=TZ))
ax.grid(alpha=0.3); ax.legend(loc='lower right',fontsize=8,ncol=2)
fig.autofmt_xdate(rotation=30,ha='right'); fig.tight_layout()
fig.savefig('/tmp/score_progress.png',dpi=150)
print('points',len(merged),'runs',len(series))
for k,v in series.items(): print(f'{k}: {len(v)} 点, 最后 {v[-1][1]:,} @ {v[-1][0]:%m-%d %H:%M}')
print('envelope end',env[-1])
