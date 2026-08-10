# -*- coding: utf-8 -*-
"""Crossref 批量核验 v54 参考文献 44 条 + 提取完整作者（AJCR 禁 et al.）"""
import re, json, time, urllib.request, urllib.parse

PROXY = 'http://127.0.0.1:7897'
handler = urllib.request.ProxyHandler({'http': PROXY, 'https': PROXY})
opener = urllib.request.build_opener(handler)
UA = {'User-Agent': 'ref-verify/1.0 (mailto:research@gmc.edu.cn)'}

SRC = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0/v54_modified_text.txt'
OUT = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0/审查与修订记录/v54_参考文献核验_2026-08-08.json'

# 提取 References 条目
refs = []
with open(SRC, encoding='utf-8') as f:
    for line in f:
        m = re.match(r'^\[(?:1[0-9]{2}|3[6-9]|4[0-4])\]\s+(\d+)\.\s+(.*)$', line.strip())
        if m:
            refs.append({'num': int(m.group(1)), 'text': m.group(2).strip()})
refs.sort(key=lambda x: x['num'])
print('提取参考文献条目:', len(refs))

def get_doi(text):
    m = re.search(r'doi:\s*([0-9]{2}\.[0-9]+/[^\s]+)', text, re.I)
    if m:
        return m.group(1).rstrip('.')
    # 尝试 DOI URL
    m2 = re.search(r'https?://(?:dx\.)?doi\.org/([^\s]+)', text)
    return m2.group(1).rstrip('.') if m2 else None

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
                if given and family:
                    authors.append(f'{family} {given}')
                elif family:
                    authors.append(family)
            return {
                'ok': True,
                'doi': doi,
                'title': (m.get('title') or [''])[0],
                'journal': (m.get('container-title') or [''])[0],
                'year': m.get('issued', {}).get('date-parts', [[None]])[0][0],
                'volume': m.get('volume', ''),
                'page': m.get('page', m.get('article-number', '')),
                'authors': authors,
            }
        except Exception as e:
            if attempt < 2:
                time.sleep(3)
            else:
                return {'ok': False, 'doi': doi, 'error': str(e)[:120]}

results = []
for r in refs:
    doi = get_doi(r['text'])
    if not doi:
        results.append({'num': r['num'], 'ok': False, 'error': 'NO_DOI_IN_TEXT', 'text': r['text'][:150]})
        print(f'[{r["num"]}] NO DOI:', r['text'][:80])
        continue
    cr = crossref(doi)
    cr['num'] = r['num']
    cr['orig_text'] = r['text']
    results.append(cr)
    status = 'OK' if cr['ok'] else 'FAIL'
    n_auth = len(cr.get('authors', []))
    print(f"[{r['num']}] {status} | {cr.get('title', '')[:50]} | authors={n_auth}")
    time.sleep(0.5)

with open(OUT, 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=1)
print()
ok = sum(1 for r in results if r['ok'])
print(f'核验完成: {ok}/{len(results)} OK, 结果已保存 {OUT}')
