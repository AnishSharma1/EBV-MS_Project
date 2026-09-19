from pathlib import Path
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, KeepTogether
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_LEFT
import json, fitz

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/pdf'; FIG=ROOT/'outputs/manuscript_figures_only_2026-09-06'
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Body',fontName='Times-Roman',fontSize=11,leading=15,spaceAfter=9))
styles.add(ParagraphStyle(name='Head',fontName='Helvetica-Bold',fontSize=14,leading=18,textColor=HexColor('#173D5B'),spaceAfter=12))
styles.add(ParagraphStyle(name='Sub',fontName='Helvetica-Bold',fontSize=10.5,leading=14,spaceBefore=7,spaceAfter=6))
styles.add(ParagraphStyle(name='Cap',fontName='Times-Roman',fontSize=9,leading=12,spaceAfter=10))
styles.add(ParagraphStyle(name='Note',fontName='Helvetica',fontSize=9,leading=12,textColor=HexColor('#725021'),spaceAfter=12))
styles.add(ParagraphStyle(name='Title2',fontName='Times-Bold',fontSize=22,leading=26,spaceAfter=15))
story=[]; source=[]
def p(text,style='Body'):
 story.append(Paragraph(text,styles[style]));source.append(text+'\n')
def page(title):
 if story:story.append(PageBreak())
 p(title,'Head');source.append('\n')
def fig(name,caption,maxheight=360):
 im=Image(str(FIG/(name+'.png')))
 scale=min(468/im.imageWidth,maxheight/im.imageHeight)
 im.drawWidth=im.imageWidth*scale;im.drawHeight=im.imageHeight*scale
 story.append(KeepTogether([im,Spacer(1,8),Paragraph(caption,styles['Cap'])]));source.append('[Figure: '+name+']\n'+caption+'\n')

p('Computational assessment of EBV-self peptide similarities across HLA class II contexts','Title2')
p('Manuscript draft | 7 September 2026','Note')
p('Author list, affiliations and corresponding author: to be confirmed.','Cap')
p('Abstract','Head')
p('Epstein-Barr virus (EBV) infection is strongly associated with multiple sclerosis (MS), but the mechanisms connecting infection to disease remain incompletely resolved. Molecular mimicry between viral and human antigens is one proposed explanation. We examined an existing computational dataset containing 6,400 allele-specific EBV-self peptide comparisons across four HLA-DRB1 alleles. Sequence analysis, peptide-HLA binding predictions, register comparisons and AlphaFold 3 modeling provided complementary forms of computational evidence. A subsequent evidence review examined eight preselected peptide pairs, separately from the structural shortlist. Two BALF5-TALDO1 pairs, in HLA-DRB1*13:03 and HLA-DRB1*15:01, received medium-priority status; six pairs remained on hold because of evidence gaps. The two pairs ranked 13th and 14th, respectively, among 1,600 comparisons within their individual alleles. BALF5 is an EBV DNA-replication protein, whereas TALDO1 encodes transaldolase, a metabolic enzyme. Neither lead had binding consensus for both peptide partners in the reviewed dataset. Exploratory electrostatic comparisons failed the overall control benchmark, and TCR docking did not consistently recover reference geometry. These findings identify computational peptide-pair hypotheses for further investigation. They do not establish natural antigen presentation, T-cell cross-reactivity or an MS disease mechanism.')
p('Keywords: Epstein-Barr virus; multiple sclerosis; HLA class II; peptide similarity; molecular mimicry','Cap')
p('Draft status','Sub')
p('This draft summarizes the supplied analysis and existing figure snapshots. Methods are descriptive, and unresolved reporting details are listed at the end. No new experiments, candidate ranking or statistical significance tests were performed for this manuscript.','Note')

page('1. Introduction')
p('Multiple sclerosis (MS) is an immune-mediated disease of the central nervous system that can affect motor, sensory, visual and cognitive function. Damage to myelin and underlying nerve fibers disrupts neural signaling. Myelin surrounds axons and supports efficient conduction; in MS, the affected structures include the brain, spinal cord and optic nerves [1].')
p('The association between Epstein-Barr virus (EBV) infection and MS has received strong longitudinal support. Bjornevik et al. studied more than 10 million US military personnel, 955 of whom developed MS. MS risk increased approximately 32-fold after EBV infection, whereas a comparable association was not observed for cytomegalovirus. Serum neurofilament light, a marker of neuroaxonal injury, increased after EBV seroconversion [2]. This temporal pattern supports infection preceding the observed biomarker change; it does not identify molecular mimicry as the responsible mechanism or exclude every potential confounder.')
p('Lanz et al. identified cross-reactivity between EBV EBNA1 and the human central nervous system protein GlialCAM in an antibody derived from an MS patient [3]. That finding provides a specific example of viral-self antigen recognition. It does not establish recognition of other EBV-human pairs or validate a T-cell mechanism for the candidates considered here.')
p('T-cell recognition introduces a different molecular context. A T-cell receptor encounters a composite peptide-HLA surface rather than an isolated peptide sequence. HLA-dependent presentation can therefore affect whether sequence resemblance is relevant to immune recognition. For class II HLA molecules, the binding register describes the placement of the peptide core within the groove. Predicted registers can suggest which residues may be exposed, but they do not experimentally establish receptor contacts.')
p('Protein-structure prediction offers an additional way to describe these candidate surfaces. AlphaFold 3 models diverse biomolecular complexes [4], although a predicted structure alone cannot demonstrate stable cell-surface presentation or functional T-cell activation. In this study, computational evidence was used to organize EBV-self peptide hypotheses and identify gaps requiring experimental resolution. The scope was prioritization, not demonstration of molecular mimicry.')

page('2. Materials and methods')
p('2.1 Study design and source materials','Sub')
p('This study comprised computational sequence analysis, a separate structural-analysis track and a later evidence review of selected peptide pairs. The figure source snapshots contain 6,400 sequence-comparison records, with 1,600 records for each of HLA-DRB1*03:01, HLA-DRB1*08:01, HLA-DRB1*13:03 and HLA-DRB1*15:01. Source materials included IEDB epitope records and peptides derived from canonical human protein sequences. A sequence derived from a reference protein was not treated as an experimentally established epitope solely because of its source.')
p('2.2 Analytical tracks','Sub')
p('Register-aware sequence comparisons and HLA-binding predictions were considered alongside predicted peptide-HLA structures. The original sequence ranks were retained as within-allele descriptive ranks. The later eight-candidate evidence dossier was additive: it audited preselected candidates and did not apply a new composite score to the full discovery dataset. The structural shortlist and evidence-review outcomes were consequently treated as distinct outputs (Figure 1).')
fig('Figure_1_workflow','<b>Figure 1. Study organization.</b> Sequence and structural analyses formed separate tracks. The dashed connector denotes a subsequent targeted evidence review rather than a uniform filtering step applied to all 6,400 records. Two reviewed candidates received medium-priority status and six remained on hold. Control analyses did not provide additional candidate support.',300)

page('2. Materials and methods (continued)')
p('2.3 Evidence review and interpretation','Sub')
p('The eight-candidate dossier recorded source-defined binding consensus, register agreement and exact-HLA evidence for the viral and human peptide partners. Evidence availability was summarized separately from prediction agreement. A negative consensus field indicates that the source-defined consensus criterion was not met; it is not an experimental measurement of nonbinding. Likewise, absence of a supporting database record was treated as an evidence gap rather than proof that presentation cannot occur.')
p('2.4 Exploratory control analyses','Sub')
p('Electrostatic comparisons were assessed against the archived control-gate outcomes. The manuscript reports those outcomes without introducing a new acceptance rule. Cases that could not be evaluated were retained as a separate category. TCR docking was summarized using archived coordinate-comparison metrics for three calibration runs. Two corrected runs excluded their own pMHC reference templates; the earlier run was retained with its original status. These analyses were used to assess interpretive limitations, not to establish candidate cross-reactivity.')
p('2.5 Data presentation','Sub')
p('Sequence-score distributions, evidence fields and control metrics were plotted from saved source tables. Counts and within-allele ranks are descriptive. No p-values, false-discovery rates or probabilities of mimicry are inferred from these scores. The calibration runs are not biological replicates. Missing values and nonevaluable outcomes were not converted to negative biological findings.')
p('Reporting details pending author verification','Sub')
p('The final Methods section requires a reconciled source manifest, retrieval dates, software versions and a traceable account of how the eight reviewed candidates were selected. These details have not been reconstructed from figure labels. The link between the structural track and each later candidate claim also requires explicit documentation before submission.','Note')

page('3. Results: sequence comparisons')
p('The sequence dataset contained 1,600 comparisons in each of four HLA contexts. Figure 2 displays the within-allele score distributions and highlights the two BALF5-TALDO1 records selected in the later evidence review. In HLA-DRB1*13:03, the record had a source-reported score of 0.281 and rank 13 of 1,600. In HLA-DRB1*15:01, the corresponding score was 0.314 and rank 14 of 1,600.')
p('These ranks indicate relatively high sequence scores within the respective datasets. They do not represent percentage identity, statistical significance or a combined measure of structural and binding evidence. The two HLA contexts are separate comparisons, not replicate observations of a single validated interaction.')
fig('Figure_2_sequence_rankings','<b>Figure 2. Within-allele sequence-score distributions.</b> Each histogram contains 1,600 comparison records. Dashed markers identify the BALF5-TALDO1 records in the two indicated HLA contexts. Scores and ranks are reproduced from the original sequence table; structural geometry and binding percentiles were not included in these ranks.',350)

page('3. Results: candidate evidence review')
p('The subsequent dossier reviewed eight preselected candidates. Two BALF5-TALDO1 pairs received medium-priority status, six candidates remained on hold, and none received high-priority status. This outcome pertains to the reviewed subset and is distinct from the separate V3 structural shortlist. It should not be interpreted as a complete all-feature ranking of the 6,400 comparisons.')
p('Both BALF5-TALDO1 records showed source-defined register agreement for the viral and human peptide partners. Binding support was less consistent. For HLA-DRB1*13:03, the viral partner met the binding-consensus criterion but the human partner did not. For HLA-DRB1*15:01, neither partner met that criterion. The exact-HLA IEDB evidence count was zero arms for the former pair and one arm for the latter; neither had an exact-HLA compatible ligand-evidence arm in the displayed snapshot (Figure 3).')
fig('Figure_3_candidate_evidence','<b>Figure 3. Evidence available for eight reviewed candidates.</b> Columns reproduce binding-consensus and register-agreement fields, supporting evidence counts and Stage-1 outcomes. Evidence counts refer to the two peptide arms, not patients or experiments. No pair-level cross-reactivity is established by a positive field. Medium priority denotes the source dossier recommendation, not validation.',285)
p('The incomplete binding and presentation evidence limits the interpretation of sequence resemblance. In particular, the recommendation to investigate the two pairs does not imply that both partners have already been shown to bind their nominated HLA molecules.')

page('3. Results and 4. Discussion')
p('3.3 Control analyses','Sub')
p('Across the displayed electrostatic control assessments, 15 outcomes were failures, 12 were not evaluable and three were passes. These are metric-level outcomes across control cases, not 30 independent biological samples. The overall benchmark did not support using this approach to strengthen candidate claims (Figure S1).')
p('TCR placement RMSDs were 19.85, 9.57 and 51.23 angstroms for the earlier 1YMM run and corrected 1ZGL and 2WBJ runs, respectively. Corresponding pMHC geometry RMSDs were 10.67, 9.93 and 30.92 angstroms (Figure S2). The variation across calibration cases did not establish consistent recovery of reference geometry. Favorable model-confidence values therefore were not treated as evidence of functional cross-recognition.')
p('4. Discussion','Head')
p('The principal outcome is a bounded set of computational hypotheses with an explicit evidence audit. Two BALF5-TALDO1 records combined relatively high within-allele sequence ranks with agreement in the reviewed register fields. That combination motivated their medium-priority status, but substantial uncertainty remained in binding and presentation evidence. The distinction matters: resemblance can motivate a question without answering whether an immune receptor recognizes both antigens.')
p('The control results also constrain the conclusions. Electrostatic agreement was not reliable across the evaluated controls, and TCR docking did not consistently reproduce reference geometry. Retaining these negative or inconclusive results prevents exploratory outputs from being mistaken for independent confirmation of mimicry. Model confidence, geometric resemblance and biological recognition address different questions.')
p('Several limitations remain. The later dossier covers only eight preselected candidates and does not provide a common final ranking across the discovery universe. The manuscript lacks experimental peptide-HLA binding, confirmed registers, natural presentation and functional T-cell cross-reactivity. Database evidence is incomplete, and absence of a record cannot resolve biological absence. A patient-expression result is not included here because no finalized expression panel and directly supporting result table were incorporated into this manuscript snapshot.')
p('The present evidence supports further investigation of the reviewed peptide-pair hypotheses, while leaving the proposed immune mechanism unresolved. Any connection to MS pathogenesis would require evidence beyond the computational similarities reported here.')

page('References and manuscript completion notes')
refs=[
'1. National Institute of Neurological Disorders and Stroke. Multiple Sclerosis: Hope Through Research. https://www.ninds.nih.gov/sites/default/files/2025-05/multiple-sclerosis-hope-through-research.pdf (accessed 7 September 2026).',
'2. Bjornevik K, et al. Longitudinal analysis reveals high prevalence of Epstein-Barr virus associated with multiple sclerosis. Science. 2022;375:296-301. doi:10.1126/science.abj8222.',
'3. Lanz TV, et al. Clonally expanded B cells in multiple sclerosis bind EBV EBNA1 and GlialCAM. Nature. 2022;603:321-327. doi:10.1038/s41586-022-04432-7.',
'4. Abramson J, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature. 2024;630:493-500. doi:10.1038/s41586-024-07487-w.'
]
for ref in refs:p(ref,'Cap')
p('Items requiring author completion','Sub')
for t in ['Confirm authorship, affiliations, contributions, funding, acknowledgments and competing-interest declarations.','Confirm database and tool citations against the versions actually used; the reference list above covers the verified background and AlphaFold 3 source only.','Reconcile the source manifest and selection history before converting this descriptive Methods draft into a submission version.','Confirm data/code availability, repository release and any applicable ethics statement; no publication, repository deposit or exemption is claimed here.','Review all newly drafted prose and figure legends. BioRender artwork remains in the appendix pending correction and a suitable licensed export.']:
 p('- '+t)
p('Figure placement','Sub')
p('Figure 1 accompanies study design; Figures 2 and 3 accompany their corresponding results. Control figures follow as supplementary material. The two BioRender illustrations are retained only in a draft-artwork appendix and are not presented as evidence.','Cap')

page('Supplementary figures')
fig('Figure_S1_electrostatics_controls','<b>Figure S1. Archived electrostatic control outcomes.</b> Rows identify control systems and structure sources. Pass, fail and not-evaluable categories follow the saved control-gate record. The overall control benchmark was unsuccessful; individual passing cells do not validate the candidate-ranking approach.',480)

page('Supplementary figures (continued)')
fig('Figure_S2_TCR_docking_controls','<b>Figure S2. TCR docking calibration geometry.</b> Coordinate-comparison RMSDs for three archived calibration runs. The 1ZGL and 2WBJ runs excluded their own pMHC reference templates; 1YMM is the earlier run. Lower values indicate closer coordinate agreement under the archived comparison procedure. Runs are not biological replicates, and no new pass/fail cutoff is applied.',380)

page('Draft artwork appendix: workflow illustration')
p('Not for submission. This saved BioRender illustration includes incorrect allele labels (*04:01 and *07:01) and an unsupported interaction-energy label. Use Figure 1 for the current study organization. Retained here for artwork revision only.','Note')
fig('BioRender_workflow','BioRender workflow draft. Existing attribution is retained. Available export: 72 dpi.',400)

page('Draft artwork appendix: recognition illustration')
p('Not for submission. This conceptual illustration contains inconsistent peptide-position labels. Its molecular shapes are schematic and do not establish an experimental binding register or T-cell receptor interaction. Retained here for artwork revision only.','Note')
fig('BioRender_recognition','BioRender recognition draft. Existing attribution is retained. Available export: 72 dpi.',400)

def footer(c,d):
 c.setStrokeColor(HexColor('#D9E0E5'));c.line(72,47,540,47)
 c.setFont('Helvetica',8);c.setFillColor(HexColor('#586C7B'))
 c.drawString(72,34,'EBV-self peptide analysis | Working manuscript draft')
 c.drawRightString(540,34,str(d.page))
pdf=OUT/'EBV_MS_manuscript_draft.pdf'
SimpleDocTemplate(str(pdf),pagesize=(612,792),rightMargin=72,leftMargin=72,topMargin=54,bottomMargin=62,title='EBV-self peptide analysis - manuscript draft',author='Authorship pending confirmation').build(story,onFirstPage=footer,onLaterPages=footer)
(OUT/'EBV_MS_manuscript_draft.txt').write_text('\n'.join(source))
doc=fitz.open(pdf)
for i,page in enumerate(doc):
 page.get_pixmap(matrix=fitz.Matrix(1,1)).save(ROOT/'tmp/pdfs'/f'page_{i+1:02}.png')
texts=[p.get_text() for p in doc]
assert len(doc)==12, len(doc)
assert all(len(t)>100 for t in texts)
(OUT/'manuscript_QA.json').write_text(json.dumps({'pages':len(doc),'figures':7,'main_figures':3,'supplementary_figures':2,'draft_artwork':2,'new_analyses':False,'status':'Working draft with explicitly identified reporting gaps'},indent=2))
print(pdf)
