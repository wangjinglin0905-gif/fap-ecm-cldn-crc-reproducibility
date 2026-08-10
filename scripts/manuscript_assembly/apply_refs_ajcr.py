# -*- coding: utf-8 -*-
"""生成 v54 AJCR 全作者参考文献（展开 et al.，修正 DOI）+ 写入 v54 文本"""
import re, json

BASE = r'C:/Users/Spica/WorkBuddy/Claw/deliverables/FAP_SDC4_CD44_v5.0'
VERIFY = BASE + r'/审查与修订记录/v54_参考文献核验_修正版_2026-08-08.json'
SRC = BASE + '/v54_modified_text.txt'

data = json.load(open(VERIFY, encoding='utf-8'))
by_num = {d['num']: d for d in data}

def fmt_authors(auths, max_show=10):
    """AJCR 风格：完整作者列表（超长时前若干 + et al.；但 AJCR 禁 et al. 需全列——按完整列出）"""
    return ', '.join(auths)

def fmt_entry(d):
    n = d['num']
    orig = d['orig']
    dd = d.get('data', {})
    if not dd or not dd.get('ok'):
        # 保留原文（未核验的）
        return f'{n}. {orig}'
    authors = fmt_authors(dd.get('authors', []))
    title = dd.get('title', '')
    journal = dd.get('journal', '')
    year = dd.get('year', '')
    vol = dd.get('volume', '')
    page = dd.get('page', '')
    doi = dd.get('doi', '')
    # 构造：N. Author1, Author2, ... Title. Journal. Year;Vol:Pages. doi:xxx
    s = f'{n}. {authors}. {title}. {journal}'
    parts = []
    if year:
        parts.append(str(year))
    if vol:
        parts.append(vol)
    if page:
        parts.append(page)
    if parts:
        s += '. ' + ';'.join(parts)
    if doi and 'None' not in doi and doi != 'None (PMID-based)':
        s += f'. doi:{doi}'
    return s

# 读取 v54 文本，替换 References 区
lines = open(SRC, encoding='utf-8').read().split('\n')
out = []
replaced = set()
for line in lines:
    m = re.match(r'^\[([^\]]+)\]\s+(.*)$', line.strip())
    if m:
        idx, text = m.group(1), m.group(2)
        is_ref = (idx.isdigit() and 118 <= int(idx) <= 153) or (idx in ('36','37','38','39','40','41','42','43','44'))
        if is_ref:
            m2 = re.match(r'^(\d+)\.\s', text)
            if m2:
                n = int(m2.group(1))
                if n in by_num and by_num[n].get('data', {}).get('ok'):
                    new_entry = fmt_entry(by_num[n])
                    out.append(f'[{idx}] {new_entry}')
                    replaced.add(n)
                    continue
    out.append(line)

print('已展开作者条目:', len(replaced), '/ 44')
with open(SRC, 'w', encoding='utf-8') as f:
    f.write('\n'.join(out))

# 验证：无 et al. 残留（参考文献区）
content = open(SRC, encoding='utf-8').read()
n_etal = len(re.findall(r'^\s*\[\d+\] \d+\. .*et al\.', content, re.M))
print('References 区 et al. 残留:', n_etal)
# 抽查 Ref 36/8/1
for n in (1, 8, 36, 44):
    m = re.search(rf'^\[[^\]]+\] {n}\. (.*)$', content, re.M)
    if m:
        print(f'Ref {n}:', m.group(1)[:120])
