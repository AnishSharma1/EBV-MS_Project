from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.section import WD_SECTION
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.text import WD_BREAK
from pathlib import Path

OUT = Path('outputs/Olivia_Thomas_EBV_MS_Master_Brief_2026-09-04.docx')
OUT.parent.mkdir(exist_ok=True)

NAVY = '18324A'; BLUE = '2E6F95'; TEAL = '2B7A78'; PALE = 'EAF2F6'
GOLD = '9A6A00'; RED = '9B2C2C'; GRAY = '5B6570'; LIGHT = 'F4F6F8'; WHITE = 'FFFFFF'

doc = Document()
sec = doc.sections[0]
sec.top_margin = Inches(0.78); sec.bottom_margin = Inches(0.72)
sec.left_margin = Inches(0.82); sec.right_margin = Inches(0.82)
sec.header_distance = Inches(0.35); sec.footer_distance = Inches(0.35)

def font(run, size=None, bold=None, color=None, italic=None, name='Aptos'):
    run.font.name = name
    run._element.get_or_add_rPr().rFonts.set(qn('w:ascii'), name)
    run._element.get_or_add_rPr().rFonts.set(qn('w:hAnsi'), name)
    if size: run.font.size = Pt(size)
    if bold is not None: run.bold = bold
    if italic is not None: run.italic = italic
    if color: run.font.color.rgb = RGBColor.from_string(color)
    return run

styles = doc.styles
normal = styles['Normal']; normal.font.name = 'Aptos'; normal.font.size = Pt(10.3)
normal.paragraph_format.space_after = Pt(5); normal.paragraph_format.line_spacing = 1.13
for sty, size, color, before, after in [
    ('Title', 27, NAVY, 0, 8), ('Subtitle', 12, GRAY, 0, 14),
    ('Heading 1', 17, NAVY, 16, 7), ('Heading 2', 13.5, BLUE, 11, 5),
    ('Heading 3', 11.5, TEAL, 8, 3)]:
    s = styles[sty]; s.font.name = 'Aptos Display' if sty != 'Normal' else 'Aptos'
    s.font.size = Pt(size); s.font.color.rgb = RGBColor.from_string(color)
    s.font.bold = sty != 'Subtitle'; s.paragraph_format.space_before = Pt(before)
    s.paragraph_format.space_after = Pt(after); s.paragraph_format.keep_with_next = True

for name, color, fill in [('Key Takeaway', NAVY, PALE), ('Caution', RED, 'FCECEC')]:
    if name not in styles:
        s = styles.add_style(name, WD_STYLE_TYPE.PARAGRAPH)
    else: s = styles[name]
    s.font.name = 'Aptos'; s.font.size = Pt(10.3); s.font.color.rgb = RGBColor.from_string(color)
    s.paragraph_format.left_indent = Inches(.18); s.paragraph_format.right_indent = Inches(.12)
    s.paragraph_format.space_before = Pt(5); s.paragraph_format.space_after = Pt(7)
    s.paragraph_format.line_spacing = 1.1

def shade(cell, fill):
    tcPr = cell._tc.get_or_add_tcPr(); shd = tcPr.find(qn('w:shd'))
    if shd is None: shd = OxmlElement('w:shd'); tcPr.append(shd)
    shd.set(qn('w:fill'), fill)

def cell_margins(cell, top=90, start=110, bottom=90, end=110):
    tcPr = cell._tc.get_or_add_tcPr(); mar = tcPr.first_child_found_in('w:tcMar')
    if mar is None: mar = OxmlElement('w:tcMar'); tcPr.append(mar)
    for m, v in [('top', top), ('start', start), ('bottom', bottom), ('end', end)]:
        node = mar.find(qn('w:'+m))
        if node is None: node = OxmlElement('w:'+m); mar.append(node)
        node.set(qn('w:w'), str(v)); node.set(qn('w:type'), 'dxa')

def set_cell_text(cell, text, bold=False, color=None, size=9.1):
    cell.text = ''
    p = cell.paragraphs[0]; p.paragraph_format.space_after = Pt(0); p.paragraph_format.line_spacing = 1.05
    font(p.add_run(str(text)), size=size, bold=bold, color=color)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER; cell_margins(cell)

def table(headers, rows, widths=None, header_fill=NAVY):
    t = doc.add_table(rows=1, cols=len(headers)); t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False; t.style = 'Table Grid'
    for i, h in enumerate(headers):
        set_cell_text(t.rows[0].cells[i], h, True, WHITE, 8.8); shade(t.rows[0].cells[i], header_fill)
    for ridx, row in enumerate(rows):
        cells = t.add_row().cells
        for i, val in enumerate(row):
            set_cell_text(cells[i], val, False, NAVY if i == 0 else None, 8.7)
            if ridx % 2: shade(cells[i], LIGHT)
    if widths:
        for row in t.rows:
            for c, w in zip(row.cells, widths): c.width = Inches(w)
    doc.add_paragraph().paragraph_format.space_after = Pt(1)
    return t

def bullet(text, level=0):
    p = doc.add_paragraph(style='List Bullet' if level == 0 else 'List Bullet 2')
    p.paragraph_format.space_after = Pt(3); p.paragraph_format.left_indent = Inches(.28 + .22*level)
    p.paragraph_format.first_line_indent = Inches(-.18)
    font(p.add_run(text), size=10.1)
    return p

def numbered(text):
    p = doc.add_paragraph(style='List Number'); p.paragraph_format.space_after = Pt(4)
    font(p.add_run(text), size=10.1); return p

def para(text='', bold_lead=None, style=None):
    p = doc.add_paragraph(style=style)
    if bold_lead and text.startswith(bold_lead):
        font(p.add_run(bold_lead), bold=True)
        font(p.add_run(text[len(bold_lead):]))
    else: font(p.add_run(text))
    return p

def callout(label, text, caution=False):
    p = doc.add_paragraph(style='Caution' if caution else 'Key Takeaway')
    shade_p = OxmlElement('w:shd'); shade_p.set(qn('w:fill'), 'FCECEC' if caution else PALE)
    p._p.get_or_add_pPr().append(shade_p)
    font(p.add_run(label + '  '), bold=True, color=RED if caution else NAVY)
    font(p.add_run(text), color=RED if caution else NAVY)
    return p

def page_break(): doc.add_page_break()

# Header/footer
hp = sec.header.paragraphs[0]; hp.alignment = WD_ALIGN_PARAGRAPH.RIGHT
font(hp.add_run('EBV–MS PROJECT  |  OLIVIA THOMAS CALL'), size=8.3, bold=True, color=GRAY)
fp = sec.footer.paragraphs[0]; fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
font(fp.add_run('Prepared 3 September 2026  •  Meeting use only'), size=8, color=GRAY)

# Cover
p = doc.add_paragraph(); p.paragraph_format.space_before = Pt(48); p.paragraph_format.space_after = Pt(8)
font(p.add_run('MEETING MASTER BRIEF'), size=10, bold=True, color=TEAL)
p = doc.add_paragraph(style='Title'); font(p.add_run('EBV–MS Molecular Mimicry Project'), size=27, bold=True, color=NAVY, name='Aptos Display')
p = doc.add_paragraph(style='Subtitle'); font(p.add_run('Recent expansion, control filling, and pMHC electrostatics — what was done, what it means, and what it does not prove'), size=12, color=GRAY)
callout('ONE-SENTENCE STATUS', 'We expanded from a narrow lead set to a broader four-allele, control-contextualized candidate program; the sequence screen retained a small experimental shortlist, while the control-first electrostatics benchmark failed its locked gate and therefore electrostatics was retired from candidate ranking.')
table(['Workstream', 'What changed', 'Current scientific status'], [
    ('Expansion', 'Broader HLA and candidate universe; independent predictor/evidence review', 'Hypothesis generation; 2 medium-priority candidates, 6 on hold in the eight-candidate dossier'),
    ('Control filling', '12 candidate panels × 25 score-blind exact-HLA N3 comparison pairs', 'Ranking context complete; definitive validation still not evaluable'),
    ('Electrostatics', 'Dense near-surface pMHC field benchmarked on development controls', 'Control gate failed; no candidate-ranking use allowed')
], [1.25, 2.55, 2.95])
doc.add_paragraph()
para('Audience: Olivia Thomas  |  Call date: 4 September 2026  |  Prepared for: Anish Sharma', style='Subtitle')

page_break()
doc.add_heading('1. The big picture in plain language', level=1)
para('The project asks whether an EBV peptide and a human nervous-system peptide could look sufficiently alike when displayed by the same HLA class II molecule to justify experimental testing for immune cross-reactivity. “Look alike” is deliberately split into separate evidence layers: sequence resemblance, predicted HLA binding/register, modeled shape, modeled electrostatics, and ultimately laboratory T-cell behavior.')
callout('CORE LOGIC', 'Each computational layer is a filter, not a verdict. The more independent layers agree, the more defensible it is to spend laboratory effort on a pair. None of the computational layers alone demonstrates molecular mimicry.')
doc.add_heading('The three recent moves', level=2)
numbered('Expansion: stop overfocusing on one allele or a tiny number of attractive pairs; build a wider, frozen hypothesis universe across HLA-DRB1*03:01, *08:01, *13:03, and *15:01.')
numbered('Control filling: compare each candidate with fair, exact-HLA reference pairs selected without seeing the candidate’s score, so a rank has context rather than being an isolated number.')
numbered('Electrostatics: ask whether the electric field above the peptide–HLA surface adds useful information beyond sequence and geometry—but first require the method to recover known development controls.')
doc.add_heading('The project’s evidence ladder', level=2)
table(['Layer', 'Question answered', 'What it cannot establish'], [
    ('Sequence', 'Do EBV and self peptide residues resemble one another?', 'HLA presentation or immune recognition'),
    ('Binding/register', 'Could both bind the same HLA, and in what P1–P9 alignment?', 'TCR binding or activation'),
    ('Geometry', 'Do modeled exposed positions occupy similar 3-D locations?', 'A real TCR footprint or induced fit'),
    ('Electrostatics', 'Do modeled near-surface electric fields resemble one another?', 'Specificity, cross-reactivity, or disease mechanism'),
    ('Experiment', 'Does the same T-cell receptor respond to both pMHCs?', 'By itself, population-level MS causation')
], [1.05, 2.55, 3.15])

doc.add_heading('2. Expansion: what was broadened', level=1)
doc.add_heading('A. Four-allele structural universe', level=2)
para('The work moved beyond a DRB1*15:01-only story. A frozen AlphaFold Server V2 universe covered four HLA-DRB1 alleles. Each job modeled one peptide–HLA complex with five AlphaFold models. Download audits treated Finder-suffixed copies as technical duplicates, not extra replicates, and checked the request identity plus all model/confidence/full-data files.')
bullet('Frozen V2 universe: 320 unique jobs across four HLA-DRB1 alleles.')
bullet('The audit trail progressed from 239 unique complete jobs, to 301/320, then 314/320 with six exact missing jobs isolated for salvage.')
bullet('The important habit: completeness was based on canonical request identity and bundle contents, not folder count.')
callout('HOW TO SAY IT', '“We expanded the structural comparison universe across four HLA risk contexts and maintained a strict manifest so duplicate downloads could not masquerade as independent models.”')

doc.add_heading('B. High-yield candidate expansion', level=2)
para('Ten fresh HLA-specific hypotheses were selected additively from the existing V3 universes. The old rankings were not changed. Eligibility required both original binding percentiles ≤20, complete five-model ensembles, and a uniquely contained declared register. Diversity rules prevented repeatedly choosing nested cores or the same EBV/self protein families.')
table(['Item', 'Count', 'Outcome'], [
    ('Fresh candidates', '10', 'All held for register or binding conflict after independent screening'),
    ('Existing BALF5–TALDO1 lead rechecks', '2', 'Both advanced with caution'),
    ('Peptide arms', '24', 'All predictor records complete'),
    ('Independent predictor records', '96', 'NetMHCIIpan 4.3 and MixMHC2pred 2.1 kept separate')
], [2.7, .7, 3.35])
para('“Advance with caution” meant only that both peptide arms retained independent binding support and both predictors agreed with the declared register. It did not mean presentation or cross-reactivity had been demonstrated.')

doc.add_heading('C. Eight-candidate evidence dossier', level=2)
para('The eight sequence-supported candidates were then audited across independent binding predictors, IEDB assay provenance, immunopeptidome searches, human/EBV sequence rarity, conservation, and common missense variation. The stage-one assay-prioritization gate completed with 0 high-priority, 2 medium-priority, and 6 hold candidates. Stage two remained not evaluable pending experimental binding and register results.')
callout('WHY THIS MATTERS', 'The expansion did not manufacture more “winners.” It increased the search breadth and then allowed most candidates to be downgraded when independent evidence or register agreement was weak. That is a strength, not a disappointment.')

page_break()
doc.add_heading('3. Control filling: giving ranks a fair comparison', level=1)
para('A raw similarity score is hard to interpret. Control filling surrounded each frozen candidate with a comparison panel: 25 exact-HLA N3 pairs chosen without using geometry or electrostatics. The candidate plus 25 comparators created a 26-pair panel, so a rank of 1/26 means the candidate looked most similar under that metric within its panel.')
doc.add_heading('What was completed', level=2)
table(['Quantity', 'Verified value'], [
    ('Frozen candidate panels', '12/12 complete'),
    ('N3 comparison pairs per panel', '25'),
    ('Total panel rows', '312 = 12 targets + 300 comparisons'),
    ('Targets ranking in top three on the V3 primary panel metric', '8/12'),
    ('New untouched strict positive-control systems admitted', '0'),
    ('Discovery unlock', 'False')
], [3.8, 2.95])
doc.add_heading('What N3 means', level=2)
callout('N3 ≠ NEGATIVE', 'N3 means recognition status is unknown. These pairs are computational comparators that help define a within-HLA background distribution. They are not confirmed non-cross-reactive biological negatives.')
para('The panels were score-blind: eligibility and ordering used peptide length, binding balance, frozen identifiers, and a deterministic hash—not the similarity outcome. This reduces cherry-picking. The originally proposed global exclusion of every target arm was infeasible because too few unique self arms remained, so the locked resolution excluded the current target’s arms and all confirmed-control ligands while flagging any comparator that reused an arm from another high-yield target.')

doc.add_heading('Three control categories you must not mix up', level=2)
table(['Category', 'Role in this project', 'Available?'], [
    ('N3 comparison pair', 'Defines the computational background for ranking; recognition unknown', 'Yes: 25 per candidate panel'),
    ('Development positive control', 'Known literature system used while developing/testing a method', 'Yes: Hy.2E11, Ob.1A12, Hy.1B11'),
    ('Untouched external validation control', 'Independent known system never used in method development', 'No adequate set; 0 newly admissible systems')
], [1.65, 3.7, 1.4])
para('Because there were zero newly admissible untouched strict systems, definitive validation was not evaluable. The three known systems remained development controls. The specificity gate was also not evaluable because no explicit N1/N2 negative registry existed.')
callout('THE HONEST CONCLUSION', 'The filled panels provide useful HLA-specific rank context. They do not validate the scoring system biologically, establish specificity, or unlock discovery claims.')

doc.add_heading('4. Why register handling dominates the interpretation', level=1)
para('An HLA-II peptide can slide in the groove. The nine-residue P1–P9 core determines which residues point into HLA pockets and which face upward toward a potential TCR. Two peptides can appear similar in sequence but present different upward-facing patterns if their registers differ.')
bullet('Register lane: declared binding cores were robust enough to compare as the primary structural interpretation.')
bullet('Sequence lane: sequence-supported candidates whose declared cores were not robust; geometry and electrostatics were retained only as sensitivity diagnostics.')
para('In the 12 filled panels, all eight top-three target results came from sequence-lane candidates with unresolved registers, while the four register-robust candidates did not place in the top three. That is why an exciting rank could not be promoted into formal support.')
callout('MEMORIZE', '“A strong rank with an uncertain register is a clue, not a result. If the peptide is shifted in the groove, the TCR-facing surface we modeled may be the wrong one.”')

doc.add_heading('5. Electrostatics: what we measured', level=1)
para('Electrostatics describes how positive and negative charges shape the electric field near the solvent-exposed pMHC surface. A TCR approaches this charged surface, so electrostatic resemblance could in principle supplement sequence and shape. But the calculation does not include a real TCR, full molecular dynamics, induced fit, detailed water networks, or experimental binding.')
doc.add_heading('Candidate pilot and six-candidate expansion', level=2)
para('The pilot tested two BALF5–TALDO1 candidates. A later package tested six additional sequence-supported candidates. Each target was ranked within its frozen exact-HLA 26-pair panel. The primary frozen rule required a full-pMHC rank of 1–3, stability across protein dielectric values 2/4/8, and a robust register.')
table(['Candidate', 'HLA', 'Full-pMHC rank (ε=2)', 'ε=4 / 8', 'HLA-subtracted rank'], [
    ('HY03_SEQ_01', 'DRB1*03:01', '19/26', '18/26 / 15/26', '18/26'),
    ('HY03_SEQ_02', 'DRB1*03:01', '24/26', '25/26 / 25/26', '13/26'),
    ('HY08_SEQ_01', 'DRB1*08:01', '13/26', '11/26 / 10/26', '24/26'),
    ('HY08_SEQ_02', 'DRB1*08:01', '24/26', '24/26 / 24/26', '17/26'),
    ('HY13_SEQ_01', 'DRB1*13:03', '12/26', '12/26 / 12/26', '9/26'),
    ('HY13_SEQ_02', 'DRB1*13:03', '21/26', '21/26 / 21/26', '25/26'),
    ('HY15_SEQ_01', 'DRB1*15:01', '8/26', '8/26 / 8/26', '23/26'),
    ('HY15_SEQ_02', 'DRB1*15:01', '19/26', '19/26 / 18/26', '5/26')
], [1.35, 1.28, 1.35, 1.35, 1.45])
callout('CANDIDATE RESULT', 'None of the eight candidates met formal sequence-plus-electrostatics support. They remained sequence-supported hypotheses; electrostatics was not supportive within these panels. This does not prove nonrecognition.')

doc.add_heading('QC correction in the first electrostatics analysis', level=2)
para('The first target-surface sampling scheme failed a common-solvent-accessibility check, meaning some points were not comparable across all structures. That initial analysis was discarded. The corrected method used 25 position-matched field points around exposed P2/P3/P5/P7/P8 positions that passed accessibility in all 60 models per panel. The correction used geometry only and did not inspect electrostatic scores.')
doc.add_heading('Core physical setup', level=2)
bullet('Structures converted with PDB2PQR 3.7.1; PROPKA protonation at pH 7.4; PARSE parameters.')
bullet('APBS 3.4.1 solved the Poisson–Boltzmann equation on a shared panel grid.')
bullet('Primary V2 control setting: nonlinear PB, protein dielectric 4, solvent dielectric 78.5, 0.15 M monovalent salt, maximum grid spacing 0.5 Å.')
bullet('Physical and sampling sensitivities were reported separately; no weighted composite score was created.')

doc.add_heading('6. The control-first electrostatics V2 benchmark', level=1)
para('Because the candidate result alone could reflect a weak endpoint rather than truly unusual biology, the method was rebuilt as a control-first test. The earlier sparse descriptor was replaced by a dense, model-specific near-surface pMHC map. Only frozen development controls were evaluated; candidate files were not read.')
table(['Scale / check', 'Result'], [
    ('Control structures', '370 total: 360 AlphaFold + 10 experimental structures'),
    ('APBS calculations', '2,220 complete'),
    ('Locked gate requirements', '30'),
    ('Package assertions', '18/18 passed'),
    ('Analyze-stage deterministic rebuild', '12,982 files byte-identical across two rebuilds'),
    ('Candidate-path contamination', 'Absent'),
    ('Final control gate', 'FAIL — candidate evaluation prohibited')
], [3.0, 3.75])

doc.add_heading('What the development controls did', level=2)
table(['Control / layer', 'Peptide electrostatics', 'Composite electrostatics', 'Surface shape'], [
    ('Hy.2E11 BALF5–MBP / PDB', '5 (fail)', '6 (fail)', '4 (fail)'),
    ('Ob.1A12 ENGA–MBP / PDB', '6 (fail)', '6 (fail)', '5 (fail)'),
    ('Hy.1B11 UL15–MBP / AF seed 271828', '12 (fail)', '10 (fail)', '10 (fail)'),
    ('Hy.1B11 UL15–MBP / AF seed 314159', '3 (unstable)', '2 (unstable)', '6 (fail)'),
    ('Hy.2E11 BALF5–MBP / AF seed 314159', '18 (fail)', '6 (fail)', '1 (pass)'),
    ('Ob.1A12 ENGA–MBP / AF seed 314159', '1 (pass)', '1 (pass)', '10 (fail)')
], [2.8, 1.3, 1.45, 1.25])
para('Some additional rows were not evaluable because at least one member of the comparison panel failed the locked ≥90% pairwise surface-map coverage rule. A good raw rank did not override missing coverage, instability, or sensitivity failure.')
page_break()
callout('MOST IMPORTANT OUTCOME', 'The computation succeeded, but the scientific gate failed. The endpoint did not reliably recover the known development controls across required structures, seeds, sensitivities, and resampling. Therefore candidate scoring stayed locked and electrostatics was permanently retired from candidate ranking.')

doc.add_heading('7. How to interpret a “failed gate” correctly', level=1)
table(['Wrong interpretation', 'Correct interpretation'], [
    ('“The electrostatics code failed.”', 'No. The package verified successfully: 370 structures, 2,220 calculations, 18 assertions, and deterministic rebuilding.'),
    ('“The candidates are disproven.”', 'No. The endpoint failed to earn trust; it cannot be used to reject or promote candidates.'),
    ('“Electrostatics never matters biologically.”', 'Not established. This particular modeled endpoint and control design were insufficiently reliable for ranking.'),
    ('“A rank-1 control proves the method.”', 'No. The locked gate required consistent recovery across controls, models, sensitivity settings, coverage, and resampling.'),
    ('“We should tune the method until it passes.”', 'That would risk overfitting. The locked failure is retained; any future redesign must be a new, separately preregistered method.')
], [2.45, 4.3])

doc.add_heading('8. What is established vs. still unknown', level=1)
table(['Status', 'What we can say'], [
    ('Established computationally', 'A reproducible four-allele structural library and score-blind panel framework were built; sequence-supported candidates were prioritized; dense electrostatics was tested under a locked control gate.'),
    ('Negative/limiting result', 'The V2 electrostatics endpoint failed the control gate and is not allowed in candidate ranking. None of eight sequence-supported candidates gained formal electrostatic support.'),
    ('Still unknown', 'Natural presentation, the correct in-cell register, TCR binding, activation, specificity, cross-reactivity, and causal relevance to MS.'),
    ('Highest-value next evidence', 'Experimental HLA-II binding and register confirmation for the two medium-priority candidates, followed—only if supported—by T-cell assays with appropriate positive and true negative controls.')
], [1.7, 5.05])

doc.add_heading('9. Your 90-second explanation to Olivia', level=1)
para('“Recently I expanded the project from a narrow set of leads into a four-allele, register-aware candidate framework. I kept the old discovery rankings frozen and used independent binding predictors and provenance checks rather than blending everything into one score. That left two medium-priority candidates for binding and register experiments and placed six on hold.”')
para('“I then filled exact-HLA comparison panels around twelve candidates. Each candidate was ranked against twenty-five score-blind N3 pairs. Eight candidates ranked in the top three on the earlier computational metric, but those were sequence-lane cases with unresolved registers. Also, N3 means recognition unknown—not biological negative—so this provides rank context, not specificity validation.”')
para('“Finally, I tested whether surface electrostatics added useful information. The candidate pilot was not supportive, so I built a control-first dense electrostatics benchmark using three known development-control systems. The computational package completed and was reproducible, but the locked control gate failed because known controls were not recovered consistently across structures, seeds, sensitivities, coverage, and resampling. Therefore I retired electrostatics from candidate ranking rather than tuning it after seeing the result. The next decisive step is experimental peptide–HLA binding and register confirmation, then T-cell testing if those pass.”')

doc.add_heading('10. Likely questions and strong answers', level=1)
qa = [
('Why these four HLA alleles?', 'They broaden the study beyond DRB1*15:01 while retaining exact-allele comparisons. Alleles are analyzed separately because the HLA groove changes both binding register and the presented surface.'),
('Why 25 comparison pairs?', 'It creates a manageable, frozen within-HLA reference panel and yields an interpretable 1–26 rank. It is still a descriptive panel, not a population-calibrated false-discovery estimate.'),
('Why call them N3?', 'Their recognition status is unknown. The label prevents us from accidentally treating untested pairs as confirmed biological negatives.'),
('Why are the strong ranks not enough?', 'Most occurred in register-uncertain sequence lanes. If the register is wrong, the modeled TCR-facing residues and surface comparison may also be wrong.'),
('Why abandon electrostatics?', 'We did not abandon electrostatics as a biological concept; we retired this endpoint from ranking because it failed the predeclared control benchmark.'),
('What would rescue the project?', 'The main project does not need rescue. Electrostatics was supplemental. Binding/register experiments can test the sequence-prioritized hypotheses directly, after which T-cell assays become interpretable.'),
('What is the strongest methodological choice?', 'Freezing rules before reading outcomes, using score-blind comparators, retaining failed gates, and keeping sequence, geometry, electrostatics, and immune-function claims separate.')]
for q, a in qa:
    p = doc.add_paragraph(); p.paragraph_format.space_before = Pt(4); p.paragraph_format.space_after = Pt(2)
    font(p.add_run('Q: ' + q), bold=True, color=NAVY)
    p = doc.add_paragraph(); p.paragraph_format.left_indent = Inches(.22); p.paragraph_format.space_after = Pt(5)
    font(p.add_run('A: ' + a))

page_break()
doc.add_heading('11. Questions worth asking Olivia', level=1)
for x in [
    'Is the distinction between register-robust and sequence-only lanes biologically sensible, or should uncertain-register candidates be excluded earlier?',
    'Are N3 exact-HLA pairs a defensible background for descriptive ranking, and what would she consider an adequate true-negative control set?',
    'For the two medium-priority pairs, what is the minimum experimental sequence: competitive HLA binding, register mapping, tetramers, activation assays, or another order?',
    'Should BALF5–MBP, ENGA–MBP, and UL15/PMM–MBP remain development controls, or are there better independent HLA-II cross-reactivity systems?',
    'Does retiring electrostatics from ranking strengthen the study’s credibility, and is there any narrower descriptive role worth retaining?',
    'Which claim boundary would she use in a poster: “computational prioritization,” “structural resemblance,” or something even narrower?'
]: bullet(x)

doc.add_heading('12. Active-recall check: answer before the call', level=1)
checks = [
('What is the single biggest difference between a sequence-lane and a register-lane candidate?', 'Register-lane candidates have a robust declared P1–P9 binding register; sequence-lane candidates do not, so surface results are sensitivity evidence only.'),
('Why is N3 not a negative control?', 'Because TCR recognition was never established as absent; it is unknown.'),
('What were the control-filling totals?', '12 panels, 25 N3 pairs each, 312 rows total, 8/12 top-three target ranks.'),
('What happened to the ten fresh expansion candidates?', 'All ten were held for register or binding conflicts after independent screening; the two existing BALF5–TALDO1 rechecks advanced with caution.'),
('What was the electrostatics candidate result?', '0/8 gained formal sequence-plus-electrostatics support.'),
('What was the electrostatics control result?', 'The package verified computationally, but the locked scientific gate failed; candidate evaluation was prohibited and electrostatics was retired from ranking.'),
('What can the project still do next?', 'Run experimental HLA binding and register confirmation on the two medium-priority pairs, then proceed to T-cell assays only if supported.')]
for i, (q, a) in enumerate(checks, 1):
    p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(1)
    font(p.add_run(f'{i}. {q}'), bold=True, color=NAVY)
    p = doc.add_paragraph(); p.paragraph_format.left_indent = Inches(.24); p.paragraph_format.space_after = Pt(5)
    font(p.add_run('Answer: '), bold=True, color=TEAL); font(p.add_run(a))

doc.add_heading('13. Final claim boundary', level=1)
callout('SAFE LANGUAGE', 'This work prioritizes peptide–HLA hypotheses and evaluates model-derived sequence, geometry, and surface-field resemblance. It does not establish natural presentation, TCR recognition or binding, activation, specificity, cross-reactivity, molecular mimicry, probability, false-discovery rate, or an MS mechanism.', caution=True)

page_break()
doc.add_heading('Appendix A. Verified artifact trail', level=1)
refs = [
('Expansion', 'processed/high_yield_candidate_expansion_2026-08-28'),
('Evidence dossier', 'processed/high_yield_candidate_evidence_2026-08-28'),
('Control-filled ranking context', 'processed/high_yield_control_validation_2026-08-28'),
('Electrostatics pilot', 'processed/pmhc_surface_electrostatics_pilot_2026-08-29'),
('Six-candidate electrostatics expansion', 'processed/pmhc_surface_electrostatics_sequence_expansion_2026-08-30'),
('Eight-candidate joined electrostatics summary', 'processed/pmhc_surface_electrostatics_all_sequence_candidates_2026-08-30'),
('Control-first electrostatics V2', 'processed/pmhc_surface_electrostatics_v2_controls_2026-08-30')]
table(['Workstream', 'Authoritative project-relative path'], refs, [2.25, 4.5])
para('Authoritative project root: /Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication', style='Subtitle')
para('The control-first electrostatics code and compact reproducibility package were also published to the Science-Fair-Project repository at commit 0de06976f4da092edf49af8898cf4420b6a678a0. Multi-gigabyte raw solver intermediates remain in the authoritative iCloud archive.')

# Avoid orphaned headings and add page-number field in footer.
for p in doc.paragraphs:
    if p.style and p.style.name.startswith('Heading'):
        p.paragraph_format.keep_with_next = True

doc.save(OUT)
print(OUT.resolve())
