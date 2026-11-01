# WP4 pipeline audit, step 2 of 3 (amendment A3): appends the round-2 extended-search records (47) and studies A11-A12.
# Writes audit_screening_log.csv and audit_table.csv in the current folder. Run order:
# wp4_make_audit_tables.py -> wp4_a3_update.py -> wp4_spotcheck_update.py
import csv
Q = {1:'"Healthy Brain Network" EEG preprocessing', 2:'"HBN-EEG"', 3:'"HBN" EEG resting state',
     4:'ds005505 OR ... OR ds005516 EEG', 5:'"EEG Foundation Challenge" HBN'}
D='duplicate'; E='exclude'
r2 = [
(1,1,'arXiv 2203.02552',D,'round 1'),(1,2,'GitHub nemarDatasets/on005507',E,'dataset page'),
(1,3,'arXiv 2503.19924',D,'round 1'),(1,4,'arXiv 2011.12107 (EEG-GCNN)',E,'not HBN (TUAB, LEMON per record)'),
(1,5,'PMC11994696',D,'round 1'),(1,6,'arXiv 2506.22455',D,'round 1'),(1,7,'arXiv 2511.15393',D,'round 1'),
(1,8,'arXiv 2606.30104',D,'round 1'),
(1,9,'ResearchGate 355928130 (NeurML)','include-pending','HBN preprocessing/ML study per record; ResearchGate blocks Claude; one open-copy lookup search found none -> manual extraction'),
(2,1,'bioRxiv 615261 (pdf)',D,'round 1 (descriptor)'),(2,2,'bioRxiv 615261v1',E,'dataset descriptor'),
(2,3,'GitHub nemarDatasets/on005516',E,'dataset page'),(2,4,'GitHub nemarDatasets/on005510',E,'dataset page'),
(2,5,'OpenNeuro ds005507',E,'dataset page'),(2,6,'nemar-dev ds005505',E,'dataset page'),
(2,7,'OpenNeuro ds005515',E,'dataset page'),(2,8,'OpenNeuro ds005505 download',E,'dataset page'),
(2,9,'OpenNeuro ds005506',E,'dataset page'),
(2,10,'arXiv 2509.00314 (CoMET)','include','HBN-EEG pretraining analysis; methods accessible'),
(3,1,'bioRxiv 615261v2.full',E,'dataset descriptor'),(3,2,'ResearchGate HBN-EEG descriptor',E,'dataset descriptor'),
(3,3,'NEMAR on005515',E,'dataset page'),(3,4,'Nature Sci Data 2017 sdata2017181',E,'dataset descriptor (HBN Biobank)'),
(3,5,'NEMAR on005516',E,'dataset page'),(3,6,'OpenNeuro ds005507',E,'dataset page'),
(3,7,'PMC11691402',E,'not HBN (older adults abstract)'),(3,8,'PMC3847481',E,'review; not HBN'),
(3,9,'arXiv 2510.22364',E,'not HBN (creativity resting EEG, per record; full text not checked)'),
(3,10,'neuromechanist.github.io HBN page',E,'dataset web page'),
(4,1,'NEMAR ds005505',E,'dataset page'),(4,2,'nemar-dev ds005505',E,'dataset page'),
(4,3,'NEMAR ds005506',E,'dataset page'),(4,4,'NEMAR ds005508',E,'dataset page'),
(4,5,'NEMAR ds005507',E,'dataset page'),(4,6,'HuggingFace xw1216 HBN archives',E,'dataset archive'),
(4,7,'EEGDash OpenNeuro list',E,'dataset listing'),(4,8,'OpenNeuro ds005507',E,'dataset page'),
(4,9,'EEGDash datasets API',E,'dataset listing'),
(5,1,'arXiv 2506.19141 (html)',D,'round 1 (challenge descriptor)'),(5,2,'ResearchGate challenge paper',E,'challenge descriptor'),
(5,3,'themoonlight.io review',E,'third-party summary'),(5,4,'fcon_1000 OpenScience',D,'round 1'),
(5,5,'vscentrum.be blog',D,'round 1'),(5,6,'MNE discourse',D,'round 1'),
(5,7,'eeg2025 proposal pdf',E,'challenge descriptor'),(5,8,'GitHub STEEGFormer',E,'code repository, not a paper'),
(5,9,'eeg2025.github.io',E,'challenge web page'),
]
with open('audit_screening_log.csv','a',newline='') as f:
    w=csv.writer(f)
    for q,r,s,d,why in r2: w.writerow([f'A3-{q}',Q[q]+' [extended]',r,s,d,why])
rows=[
('A11','arXiv 2509.00314','Li et al., 2025 (CoMET)','HBN-EEG (Shirazi et al. 2024), 3155 subjects','none','-','none described (Cz-referenced; 62 shared channels selected)','No channel-level detection','Experiment: Pre-training; Appendix 1 and 3',''),
('A12','ResearchGate 355928130','NeurML (authors not retrieved)','HBN (per record)','PENDING','PENDING','PENDING','PENDING - manual extraction','-','ResearchGate blocked; no open copy found'),
]
with open('audit_table.csv','a',newline='') as f:
    w=csv.writer(f); [w.writerow(r) for r in rows]
print(len(r2),'round-2 records;',sum(d.startswith('include') for *_,d,_ in [(a,b,c,d,e) for a,b,c,d,e in r2]),'new eligible')
