# WP4 pipeline audit, step 1 of 3: round-1 screening log (50 records) and audit table (studies A01-A10), as screened on 26 Oct 2026.
# Writes audit_screening_log.csv and audit_table.csv in the current folder. Run order:
# wp4_make_audit_tables.py -> wp4_a3_update.py -> wp4_spotcheck_update.py
import csv
Q = {1:'"Healthy Brain Network" EEG preprocessing', 2:'"HBN-EEG"', 3:'"HBN" EEG resting state',
     4:'ds005505 OR ds005506 OR ... OR ds005516 EEG', 5:'"EEG Foundation Challenge" HBN'}
# (query, rank, source/id, decision, reason)
log = [
(1,1,'arXiv 2506.22455','include','HBN-EEG resting analysis; methods accessible'),
(1,2,'arXiv 2203.02552','include','HBN resting analysis; methods accessible'),
(1,3,'arXiv 2310.11994','include','HBN resting analysis; methods accessible'),
(1,4,'PMC13655006 (Front Neuroinform 2026)','exclude','not HBN data (sleep recordings, synthetic data)'),
(1,5,'arXiv 2503.19924','include','HBN resting analysis (inattentive ADHD vs controls); methods accessible'),
(1,6,'arXiv 2502.08686','exclude','not HBN data (LEMON, semi-simulated EOG dataset); full text checked'),
(1,7,'arXiv 2511.15393','exclude','not HBN data (single-hospital cohort); full text checked'),
(1,8,'arXiv 2606.30104','include','HBN-EEG pretraining analysis; methods accessible'),
(1,9,'arXiv 2606.23706','include','HBN (challenge distribution) analysis; methods accessible'),
(1,10,'conf.msu.ru Lomonosov-2022 abstract','exclude','not HBN; conference abstract on infant EEG'),
(2,1,'bioRxiv 10.1101/2024.10.03.615261 (pdf)','exclude','dataset descriptor'),
(2,2,'arXiv 2009.13402','exclude','review'),
(2,3,'PMC10659091 (Front Hum Neurosci 2023)','include','HBN resting analysis; methods accessible'),
(2,4,'bioRxiv 615261v2 (page)','exclude','duplicate of Q2 rank 1 (descriptor)'),
(2,5,'arXiv 2310.11994','duplicate','already screened (Q1 rank 3)'),
(2,6,'arXiv 2606.30104','duplicate','already screened (Q1 rank 8)'),
(2,7,'Zenodo record 17306881','exclude','dataset record'),
(2,8,'EEGLAB newsletter issue 19','exclude','newsletter'),
(2,9,'CNS mailing-list post','exclude','project announcement'),
(2,10,'nemar-dev forum ds005505','exclude','dataset forum listing'),
(3,1,'PMC10659091','duplicate','already screened (Q2 rank 3)'),
(3,2,'arXiv 2310.11994','duplicate','already screened (Q1 rank 3)'),
(3,3,'PMC11994696 (Front Neuroimaging 2025)','include','HBN resting/video analysis; methods accessible'),
(3,4,'arXiv 2607.06629','include','HBN pretraining analysis; methods accessible'),
(3,5,'medRxiv 10.1101/2020.10.09.20207936','exclude','fMRI, not EEG'),
(3,6,'bioRxiv 10.1101/149369','exclude','dataset descriptor (HBN Biobank)'),
(3,7,'PMC12562131','include-pending','HBN resting EEG (n=100) per record; PMC full text blocked to Claude (browser check, 2 attempts) -> manual extraction needed'),
(3,8,'PMC6870341 figure page','exclude','not HBN (2012 article figure)'),
(3,9,'bioRxiv 615261v2 (page)','exclude','duplicate of Q2 rank 1 (descriptor)'),
(3,10,'ISMRM 2014 abstract 0386','exclude','not HBN'),
(4,1,'medRxiv 10.1101/2024.10.25.24316133','exclude','datasets described in record are clinical cohorts; HBN not mentioned (full text NOT checked - spot-check)'),
(4,2,'arXiv 2605.21280','exclude','not HBN (TUH corpora only); full text checked'),
(4,3,'USPTO 11529107','exclude','patent'),
(4,4,'bioRxiv 10.1101/2020.02.14.948968','exclude','not HBN; general EEG essay'),
(4,5,'USPTO 11779262','exclude','patent'),
(4,6,'USPTO 11633139','exclude','patent'),
(4,7,'alphaneural EEGDash ds005510 pointer','exclude','dataset pointer'),
(4,8,'PMC4402280 figure page','exclude','animal study'),
(4,9,'PMC5683343 table page','exclude','not HBN (Down syndrome review table)'),
(4,10,'MDPI Brain Sci 14(2):136','exclude','review; not HBN'),
(5,1,'arXiv 2506.19141 (pdf)','exclude','challenge proposal/descriptor, not an analysis'),
(5,2,'arXiv 2506.19141v1 (html)','duplicate','duplicate of Q5 rank 1'),
(5,3,'arXiv 2508.17742','exclude','not HBN (14 benchmark datasets, none HBN); full text checked'),
(5,4,'eeglablist post Aug 2025','exclude','announcement'),
(5,5,'EEGLAB newsletter issue 20','exclude','newsletter'),
(5,6,'vscentrum.be blog','exclude','blog'),
(5,7,'fcon_1000 OpenScience page','exclude','web page'),
(5,8,'MNE discourse thread','exclude','forum post'),
(5,9,'lib.upmc.fr internship posting','exclude','job posting'),
(5,10,'web3.arxiv abs 2506.19141','duplicate','duplicate of Q5 rank 1'),
]
with open('audit_screening_log.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['query_no','query','rank','record','decision','reason'])
    for q,r,s,d,why in log: w.writerow([q,Q[q],r,s,d,why])

cols=['audit_id','record','authors_year','data_distribution','channel_detection','detector_type','reference','order_category','location_in_paper','notes']
tab=[
('A01','arXiv 2506.22455','Truong & Delorme, 2025','HBN-EEG BIDS R1-R11 via EEGDash','none','-','none applied (Cz-referenced; zero Cz removed)','No channel-level detection','Methods: Preprocessing; Discussion: Data cleaning','authors state bad-channel detection (clean_rawdata) deliberately not used'),
('A02','arXiv 2203.02552','Truong et al., 2022','original HBN release (2022), resting eyes-closed','yes','clean_rawdata: flat-line 5 s, channel correlation 0.7, line noise 4; spline interpolation','linked mastoids','References before detecting','Methods: Data preprocessing','exposure limited to frozen mastoid channels; 5-s flat-line criterion alone cannot flag shorter freezes'),
('A03','arXiv 2310.11994','Hu et al., 2023','original HBN releases; Automagic','yes','Automagic recommended settings (detection, removal, interpolation)','not stated','Not reported','Methods 4.2-4.3','order relative to referencing not stated; not inferred from Automagic defaults'),
('A04','arXiv 2503.19924','Chanu et al., 2025','original HBN (diagnosis-based selection)','yes','trimOutlier amplitude range 0.0001-100 uV','no voltage re-referencing; relative phase vs global mean phase computed after rejection','Detects before referencing','Methods; Appendix B','category assigned because the only reference-like step (mean phase) follows rejection'),
('A05','arXiv 2606.30104','Yuce et al., 2026','HBN-EEG provided preprocessed version (100 Hz, 0.5-50 Hz), R1-R10','yes','noisy-channel detection + spherical spline (PREP-style citation)','average','References before detecting','Appendix A.1, Pretraining Dataset','SPOT-CHECK: sentence order could be read as average reference applied before interpolation only; most natural reading = before detection'),
('A06','arXiv 2606.23706','Jean Jacques et al., 2026','EEG2025 challenge distribution (0.5-50 Hz, 100 Hz)','none','-','none described; 128+1 channels input','No channel-level detection','Section 3.2',''),
('A07','PMC10659091','Modarres et al., 2023','HBN-provided preprocessed files (original Biobank)','yes (by HBN pipeline)','bad electrodes identified and spline-interpolated per Langer et al. 2017','not stated','Not reported','Materials and methods',''),
('A08','PMC11994696','Gomez-Tapia et al., 2025','original HBN Biobank','yes','PREP (detection, interpolation, robust average reference), then ICA','robust average (PREP)','Detects before referencing','Section 3.1; Figure 1',''),
('A09','arXiv 2607.06629','Segal et al., 2026','published derivative hbn-eeg-post-asr-ica','yes','absolute RMS > 60 uV; reference channel dropped','common average (inside ICA step)','Detects before referencing','Methods: EEG Segmentation and Preprocessing','detector is upper-threshold only; no flat-line criterion described'),
('A10','PMC12562131','not retrieved','HBN resting, n=100 (per record)','PENDING','PENDING','PENDING','PENDING - manual extraction','-','full text blocked to Claude'),
]
with open('audit_table.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(cols); [w.writerow(r) for r in tab]
print(len(log),'records logged;', sum(d.startswith('include') for _,_,_,d,_ in log),'included')
from collections import Counter
print(Counter(r[7] for r in tab))
