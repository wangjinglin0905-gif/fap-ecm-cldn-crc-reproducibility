# -*- coding: utf-8 -*-
"""v54 参考文献：修正 Freedman DOI + Kinugasa 作者 + 生成全作者列表（AJCR 禁 et al.）"""
import re, json, time, urllib.request, urllib.parse

PROXY = 'http://127.0.0.1:7897'
handler = urllib.request.ProxyHandler({'http': PROXY, 'https': PROXY})
opener = urllib.request.build_opener(handler)
UA = {'User-Agent': 'ref-verify/1.0 (mailto:research@gmc.edu.cn)'}

SRC = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0/v54_modified_text.txt'
VERIFY_JSON = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0/审查与修订记录/v54_参考文献核验_2026-08-08.json'
OUT_JSON = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0/审查与修订记录/v54_参考文献核验_修正版_2026-08-08.json'

# 1. 提取 44 条参考文献原文
refs = []
with open(SRC, encoding='utf-8') as f:
    for line in f:
        m = re.match(r'^\[([^\]]+)\]\s+(.*)$', line.strip())
        if not m:
            continue
        idx, text = m.group(1), m.group(2)
        is_ref = (idx.isdigit() and 118 <= int(idx) <= 153) or (idx in ('36','37','38','39','40','41','42','43','44'))
        if is_ref:
            m2 = re.match(r'^(\d+)\.\s(.*)$', text)
            if m2:
                refs.append({'num': int(m2.group(1)), 'orig': m2.group(2).strip()})
refs.sort(key=lambda x: x['num'])
print('参考文献:', len(refs))

# 2. 读核验结果
verified = json.load(open(VERIFY_JSON, encoding='utf-8'))
vmap = {}
for v in verified:
    if v.get('num') and v.get('ok'):
        vmap[v['num']] = v
# 修正 Ref 36 DOI
if 36 in vmap:
    vmap[36]['doi'] = '10.1080/07350015.1983.10509354'
    vmap[36]['title'] = 'A nonstochastic interpretation of reported significance levels'
    vmap[36]['journal'] = 'Journal of Business & Economic Statistics'
    vmap[36]['year'] = 1983
    vmap[36]['volume'] = '1'
    vmap[36]['page'] = '292-298'
    vmap[36]['authors'] = ['Freedman D', 'Lane D']

def crossref(doi):
    url = 'https://api.crossref.org/works/' + urllib.parse.quote(doi)
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers=UA)
            d = json.loads(opener.open(req, timeout=30).read())
            m = d['message']
            authors = []
            for a in m.get('author', []):
                given = a.get('given', '')
                family = a.get('family', '')
                if family:
                    authors.append((family + ' ' + given).strip())
            return {'ok': True, 'doi': doi, 'title': (m.get('title') or [''])[0],
                    'journal': (m.get('container-title') or [''])[0],
                    'year': m.get('issued', {}).get('date-parts', [[None]])[0][0],
                    'volume': m.get('volume', ''), 'page': m.get('page', m.get('article-number', '')),
                    'authors': authors}
        except Exception as e:
            if attempt < 2:
                time.sleep(3)
            else:
                return {'ok': False, 'doi': doi, 'error': str(e)[:100]}

# 3. Ref 8 Kinugasa（PMID 17970035）via PubMed eutils
def pubmed_authors(pmid):
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json'
    try:
        d = json.loads(opener.open(urllib.request.Request(url, headers=UA), timeout=30).read())
        res = d['result'][pmid]
        auths = [f"{a['name']}" for a in res.get('authors', [])]
        return {'ok': True, 'authors': auths, 'title': res.get('title',''),
                'journal': res.get('source',''), 'year': res.get('pubdate','')[:4],
                'volume': res.get('volume',''), 'page': res.get('pages','')}
    except Exception as e:
        return {'ok': False, 'error': str(e)[:100]}

print('查询 Ref 8 Kinugasa (PMID 17970035)...')
r8 = pubmed_authors('17970035')
if r8['ok']:
    vmap[8] = {'num': 8, 'ok': True, 'doi': 'None (PMID-based)', 'title': r8['title'],
               'journal': r8['journal'], 'year': r8['year'], 'volume': r8['volume'],
               'page': r8['page'], 'authors': r8['authors'], 'pmid': '17970035'}
    print('Ref 8 OK:', r8['authors'][:5])
else:
    print('Ref 8 FAIL:', r8['error'])

# 4. 对缺失的 Ref 36 补 Crossref（修正 DOI）
if 36 in vmap and vmap[36].get('ok'):
    print('Ref 36 使用手动修正 DOI')
elif 36 not in vmap:
    cr = crossref('10.1080/07350015.1983.10509354')
    if cr['ok']:
        cr['num'] = 36
        vmap[36] = cr
        print('Ref 36 OK via corrected DOI')

# 5. 汇总最终核验表
final = []
missing_auth = []
for r in refs:
    n = r['num']
    entry = {'num': n, 'orig': r['orig'], 'verified': n in vmap}
    if n in vmap:
        entry['data'] = vmap[n]
        if not vmap[n].get('authors'):
            missing_auth.append(n)
    final.append(entry)

print()
print('核验覆盖:', sum(1 for f in final if f['verified']), '/', len(final))
print('缺作者条目:', missing_auth or '无')

with open(OUT_JSON, 'w', encoding='utf-8') as f:
    json.dump(final, f, ensure_ascii=False, indent=1)
print('修正版核验表已保存:', OUT_JSON)
