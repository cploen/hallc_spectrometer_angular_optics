#!/usr/bin/env python3
"""Build the illustrated guide from distributed documentation figures.

Requires reportlab and Pillow. No replay data or model fitting is needed.
"""
from pathlib import Path
from xml.sax.saxutils import escape

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak

HERE=Path(__file__).resolve().parent
DOCS=HERE.parent
FIG=DOCS/'joint_surface_figures'
OUT=DOCS/'JOINT_SURFACE_METHOD_GUIDE.pdf'
WIDTH=504
NAVY=colors.HexColor('#17344c')
TEAL=colors.HexColor('#247d96')
INK=colors.HexColor('#243441')
MUTED=colors.HexColor('#52606b')
STYLES={
    'title':ParagraphStyle('title',fontName='Helvetica-Bold',fontSize=21,leading=25,textColor=NAVY,spaceAfter=15),
    'cover':ParagraphStyle('cover',fontName='Helvetica-Bold',fontSize=29,leading=34,textColor=NAVY,spaceAfter=14),
    'deck':ParagraphStyle('deck',fontName='Helvetica',fontSize=15,leading=20,textColor=TEAL,spaceAfter=14),
    'body':ParagraphStyle('body',fontName='Helvetica',fontSize=10.5,leading=14.3,textColor=INK,spaceAfter=8),
    'small':ParagraphStyle('small',fontName='Helvetica',fontSize=9.0,leading=12.0,textColor=MUTED,spaceAfter=8),
    'sub':ParagraphStyle('sub',fontName='Helvetica-Bold',fontSize=12,leading=15,textColor=TEAL,spaceBefore=8,spaceAfter=7,keepWithNext=True),
    'cell':ParagraphStyle('cell',fontName='Helvetica',fontSize=9.2,leading=12,textColor=INK),
    'head':ParagraphStyle('head',fontName='Helvetica-Bold',fontSize=9.2,leading=12,textColor=colors.white),
    'code':ParagraphStyle('code',fontName='Courier',fontSize=8.4,leading=11.5,textColor=INK,spaceAfter=8),
}
story=[]
sections=[]


def p(text,style='body'):
    story.append(Paragraph(text,STYLES[style]))


def section(title):
    if story:story.append(PageBreak())
    sections.append(title)
    p(title,'title')


def sub(title):p(title,'sub')


def fig(name,height=260,caption=None):
    path=FIG/(name+'.png')
    with PILImage.open(path) as im:w,h=im.size
    scale=min(WIDTH/w,height/h)
    obj=Image(str(path),width=w*scale,height=h*scale)
    obj.hAlign='CENTER'
    story.append(obj)
    if caption:p(caption,'small')
    story.append(Spacer(1,6))


def eq(name):fig('equation_'+name,height=42)


def table(headers,rows,widths):
    data=[[Paragraph(escape(str(v)),STYLES['head']) for v in headers]]
    data += [[Paragraph(str(v),STYLES['cell']) for v in row] for row in rows]
    t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),NAVY),('VALIGN',(0,0),(-1,-1),'TOP'),
        ('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),
        ('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f0f5f7')])]))
    story.extend([t,Spacer(1,10)])


def code(lines):p('<br/>'.join(escape(line).replace(' ','&#160;') for line in lines.splitlines()),'code')


section('Continuous band surfaces')
p('Geometric discovery, event assignment<br/>and the HMS / SHMS farm-test branch','deck')
p('A practical and mathematical guide | 27 October 2026 | Version 1')
fig('central_3d',height=282,caption='Frozen feature-only HMS pilot, run 1544, central foil, -10% to +10% delta. Measured coordinates remain unchanged. Colors are physical row/column names applied after the surface fits; gray crosses have no assignment. Each plane shows 4,800 sampled assigned events and all its rejections. This sampling emphasizes rejections and is not a density estimate.')
sub('What this version replaces')
p('The first guide described an earlier fit whose training indices largely came from reconstructed sieve positions. This guide describes the corrected method: discover ordered bands from focal-plane coordinates and continuous delta, freeze that geometry, and only then attach physical sieve names. Historical final labels are comparison material only.')
p('Branch: <b>codex/joint-surface-ifarm</b>, based on integration commit 505e656. The branch prepares conventional pre-GMM candidate trees from replay. Full-replay physics validation on ifarm remains the next experiment.','small')

section('1. The specific problem and the workflow')
p('For one foil, an X sieve row produces a broadened curve in (xfp, xpfp) at fixed delta. Across delta, the curve traces a sheet. A Y sieve column produces an analogous sheet in (yfp, ypfp, delta). Finite delta intervals superimpose changing sheet sections, which can obscure spacing, curvature and narrow gaps.')
p('The objective is reliable geometric discovery and understandable boundaries with less angle adjustment or shoulder rescue. X and Y are fitted separately. Joining an event\'s X row and Y column determines its physical sieve hole; no X-Y coupling term is required by this implementation.')
table(['Method','What establishes bands','Where physical identity enters'],[
    ['Earlier slice procedures','Angle scans, local peaks or Gaussian regions in each reporting interval. Pieces can be split deliberately.','Later alignment and merging of components.'],
    ['Preserved provisional integration','X whole-distribution curves per interval, with pooled lowest-delta geometry; Y spacing-supported peaks and Gaussian regions.','Existing relabel stages.'],
    ['This continuous method','One periodic density initialization and one ordered surface family per foil and axis, across continuous delta.','One orientation and integer origin for the frozen family.'],
],[118,205,181])
sub('What the new branch replaces')
p('It replaces the focal-plane band-discovery/relabel/candidate-construction path. It writes the usual TXCand and TYCand trees directly from event-level surface assignments. A 2D TCutG per finite delta interval cannot exactly represent a boundary that depends on each event\'s delta.')
p('Foil selection remains stage 01. Candidate output remains stage 04. GMM cleanup, the separate density cleanup, and SVD optimization are later decisions; none is silently run by the surface driver.')

section('2. Prior knowledge: exactly what is used')
table(['Stage','Prior information or assumption','Dependence on reconstruction'],[
    ['Replay admission','The existing electron PID thresholds, arm acceptance, and foil polygons in target-y versus delta.','Target-y and delta are reconstructed quantities. Incorrect foil separation or momentum reconstruction can affect the input.'],
    ['Geometry fitting','Only the axis FP position, FP slope and event delta. Smooth low-order sheets, approximate repeated spacing, and preserved ordering.','No reconstructed sieve position or old row/column label is a fit input. Delta still comes from replay.'],
    ['Hardware limit','9 possible HMS rows/columns; 11 for centered SHMS. An implausibly larger inferred span is refused.','An upper limit, not a required number of occupied bands.'],
    ['Physical naming','Known sieve layout and pitches; by default, median reconstructed sieve positions of already discovered bands choose one global map.','Can inherit a whole-family orientation/origin error. Explicit manual global maps avoid this use of reconstructed sieve positions.'],
    ['Evaluation','Old final post-GMM event pairs are read after frozen fits and names.','Comparison only. No optimization of agreement.'],
],[82,239,183])
p('The method is therefore independent of reconstructed sieve coordinates for <b>band geometry</b>. The default complete workflow is not independent of reconstruction: admission uses target-y/delta, and automatic physical naming uses sieve positions. That distinction is essential when transferring to a poorly reconstructed setting.')
p('The backup of initializing from the central row and column near delta = 0 was permitted during exploration but was not used. It is not an automatic fallback in this branch. No pre-existing event labels are needed to run it.','small')

section('3. Physical names are attached after fitting')
p('Discovery supplies an ordered integer j: 0, 1, 2, ... . It is a temporary geometric index. Zero does not automatically mean the first physical sieve row. A visible Y family with six bands can occupy six consecutive positions within a nine-column HMS grid.')
eq('classify')
p('In the naming part of this equation, s is +1 or -1 and o is one integer origin. X and Y each have their own s and o. One map applies to the entire foil/axis across delta. No individual event is moved between discovered ridges to improve agreement with the reconstruction.')
sub('Default automatic naming')
p('After the geometry files and predictions are hashed and frozen, the label stage compares each band\'s median reconstructed sieve coordinate with the nominal hardware coordinate. It chooses the valid global orientation/origin with the smallest median mismatch in units of the sieve pitch, and records alternatives. The fit coefficients and boundaries do not change.')
sub('Manual naming when the reconstruction is unreliable')
p('Supply <b>global_maps</b> in config/joint_surface.json. For example, a reviewed X family increasing with physical row, starting at physical row 1, uses orientation 1 and origin 1. All discovered accepted indices must map inside the arm\'s physical grid.')
code('"global_maps": {\n  "RUNGROUP/foil0": {\n    "x": {"orientation": 1, "origin": 1},\n    "y": {"orientation": 1, "origin": 0}\n  }\n}')
p('A manual map supplies physical names, not training labels. If a family is missing an internal position, the geometric index spacing must already represent that gap. An orientation/origin map cannot repair a missed or doubled periodicity; inspect and refit the geometry instead.')

section('4. Coordinates and a fixed reference frame')
p('For X, s = xpfp in radians and x = xfp in centimeters. For Y, substitute ypfp and yfp. Standardization subtracts each training mean and divides by its training standard deviation. Rotation produces dimensionless q and r:')
eq('rotation')
p('<b>q</b> is the standardized separation coordinate. <b>r</b> is the perpendicular rotated focal-plane coordinate. Earlier code called r "p"; that p was never momentum. The standardization and rotation are learned once per foil/axis and remain fixed across delta. The default angle is 135 degrees in standardized coordinates.')
p('Delta is stored numerically in percent: delta = -8 means an 8% offset below nominal, not -0.08. Slopes remain radians in saved event arrays; figures display milliradians. Rotating standardized coordinates does not assert a physical equality between centimeters and radians.')
fig('actual_projections',height=220,caption='Actual one-dimensional histograms from the same fixed X frame, run 1544 central foil. Each histogram includes all evaluation events in its stated interval. These are measured-coordinate projections, not unfolded model displays. Different population lengths and acceptance can change a projection without changing a physical row identity.')
p('Changing the frame independently in each delta interval could manufacture apparent tilt or curvature. This procedure keeps the frame fixed so that the model must describe the observed evolution.')

section('5. Discovering an ordered index without labels')
p('The first step seeks a smoothly varying phase coordinate whose high-density ridges occur near successive integers. A Fourier transform of the whole-foil q/r density supplies several starting spacings and directions. Numerical optimization then aligns a smooth phase across the complete event delta range.')
eq('phase')
p('Here psi is a low-order polynomial in scaled q, r and, for one candidate initializer, delta. The complex exponential repeats whenever psi changes by one. If many ridge events have the same fractional phase, their contributions add coherently and the squared magnitude C is large. The fit penalizes curvature and extra delta terms; it is not a histogram peak count forced to match nine or eleven.')
sub('From phase to provisional band indices')
p('The common fractional phase alpha is removed. Near-center events within the seed radius of an integer receive a provisional j. The observed span of these integers establishes the candidate number of positions. The surface stage then alternates coefficient fits and geometric index updates.')
sub('Why a high value is not enough')
p('Noise, a wrong periodic harmonic or an acceptance structure can also produce apparent periodicity. A separate event subset must reproduce the phase concentration. The default minimum squared phase concentration is 0.03. If neither the delta-free nor the delta-dependent initializer passes, the stage stops without claiming a discovered family.')
p('The procedure assumes approximately repeated ordered bands and low-order variation. It can miss a weak edge band or choose a harmonic. These are geometric assumptions to inspect, not truths supplied by the hardware layout.')

section('6. The family of fitted sheets')
p('Each axis uses q = F(j, r, delta). For numerical scaling define t = (j - h)/h with h = (count - 1)/2, R = r/rscale, and D = (delta - dc)/ds. rscale is the training 99th percentile of |r|. dc and ds are the midpoint and half-width of the configured delta range. Thus D = delta/10 for HMS [-10,10), and D = (delta - 6)/16 for SHMS [-10,22).')
eq('basis')
p('The displayed expression is the largest candidate. Coefficients are shared across the ordered family; there is not a separate high-order polynomial for every band. All terms are linear in their coefficients even though the surface curves in the coordinates.')
table(['Family','Terms added to the preceding family','What the added terms permit'],[
    ['bilinear','1, R, t, tR','Shared offset/slope and band-dependent spacing/tilt.'],
    ['curved','R squared, t R squared','Common and band-dependent curvature along r.'],
    ['static','t squared, t squared R, t squared R squared','Spacing and curvature that vary nonlinearly with index.'],
    ['shift','D','Shared delta-dependent displacement.'],
    ['spacing','tD','Delta-dependent separation of indices.'],
    ['tilt','RD','Delta-dependent common tilt.'],
    ['divergence','tRD','Band-dependent tilt evolution with delta.'],
],[78,181,245])
p('A static family can describe a sheet extending through delta: its visible coordinates may be strongly correlated with delta through acceptance. An explicit D term is justified only if separate model-choice events favor it. "Joint" does not mean every term must be nonzero.')

section('7. Robust fitting, regularization and ordering')
eq('objective')
p('The coefficient vector c multiplies a vector of basis terms phi. The difference q - F is a <b>focal-plane geometric residual</b>, in dimensionless q. It is not an SVD residual, a new optics-matrix residual, or a difference after replay with a refitted matrix.')
p('Repeated weighted least-squares updates reduce the influence of events far from the current centers. A diagonal coefficient penalty Lambda stabilizes poorly determined combinations. Basic terms receive a very small penalty; higher-order and explicit delta terms receive larger penalties. Population balancing prevents the densest band/delta combinations from dominating every update.')
p('Positive change with band index preserves ordering. Inequalities are checked over occupied r/delta cells and the relevant polynomial extrema. They constrain the fitted domain; they do not prove reliable ordering throughout empty r/delta space. Event inversion separately rejects locations where the fitted index direction is not increasing.')
fig('acceptance_domain',height=210,caption='Measured xfp versus delta, colored by frozen physical X index, for the complete evaluation samples in the two cases. The thin correlated domains explain why independent r and delta terms can be difficult to identify. Empty combinations are not additional measurements.')
p('Regularization can stabilize coefficients inside the observed domain. It cannot create information about an unpopulated combination of r and delta. Good prediction along the populated track does not establish a unique physical interpretation of each delta coefficient.')

section('8. Model choice and independent event roles')
p('In the exploratory study, fit A used 32% of events for geometry and 8% for model choice. Fit B used disjoint 32% and 8% subsets. A common 20% evaluation sample was used by neither fit nor model choice. The same split served both axes, and A/B shared the fixed A frame.')
p('The farm workflow uses <b>one primary fit</b>. By default, 80% of a foil\'s pre-band events determine coefficients, frame, population evidence and envelopes; 20% choose the initializer and family. The saved metadata contains the exact indices. All admitted events are then classified by the chosen frozen model. There is no additional independent physics-test subset in this default production-candidate pass.')
sub('What chooses the family')
p('The score is the negative log of a normalized conditional q density, averaged equally across populated reporting intervals. It includes ridge widths transformed from index space into q using local spacing, mixture populations and a uniform background term. It scores all model-choice events, including ones later rejected. It is not the number of returned bands, the smallest center residual, or historical-label agreement.')
p('The simplest nested family within 0.01 natural-log units per event of the best score is selected. A delta-free initializer is preferred within the same tolerance. Thus extra flexibility has to produce a measurable gain on separate events.')
sub('What the A/B comparison established')
p('A disagreement means two independently learned geometries propose different physical identities for the same evaluation event. That is different from a comparison with old final labels. A/B is a useful audit during method development; it is not required to supply two contradictory labels to SVD. The farm branch saves one explicit primary result.')

section('9. What decides a low-delta event?')
p('At the event\'s own r and delta, the model predicts q for every band shown by the separate population test to be occupied in that reporting interval. The assigned band is the closest predicted q center, provided the inversion/range and training-envelope checks also pass.')
eq('classify')
eq('boundary')
p('For two adjacent eligible centers, the dividing surface is their midpoint in q at that r and delta. This is a geometric divider, not an independently measured density valley. With a missing eligible position, nearest-center regions extend between the remaining eligible bands, subject to the other rejection rules.')
fig('low_delta_planes',height=305,caption='All evaluation events at -10% <= delta < -8%, fit A, minimum count 10. Top: central foil. Bottom: difficult -8 cm foil in runs 6983/6984. Colors are postfit physical names; gray crosses lack an assignment on that axis. These are raw measured FP planes, with no unfolding.')
p('The whole-delta fit supplies the functional shape in this low-delta section. It can resolve an overlap produced by projecting away delta or another focal-plane coordinate. It cannot identify which of two physically different populations produced an identical measured (q,r,delta). Local sparse data and independent-fit variation remain relevant.')

section('10. Population evidence is separate from geometry')
p('A fitted family predicts possible ordered positions. It does not establish that every position contains a visible population in every reporting interval. The population test uses training events in that interval, transformed to the continuous index coordinate u defined by q = F(u,r,delta).')
p('A constrained one-dimensional mixture places candidate peaks at geometric indices, fits a shared width and population weights, and includes uniform background. For each index it compares a fit including that peak with a fit omitting it, using separate events within this test. It also requires an effective population and evidence distributed along r.')
table(['Default setting','Role'],[
    ['occupancy_min_events = 100','Minimum finite training events for a reporting interval; otherwise no occupied positions are published there.'],
    ['occupancy_cap = 2000','Deterministic cap on events used for the population test.'],
    ['occupancy_likelihood_gain = 10','Required improvement of the penalized comparison including one position.'],
    ['occupancy_min_effective_events = 8','Minimum total fractional membership attributed to a position.'],
    ['occupancy_min_coherent_intervals = 2','Require near-center evidence in at least two longitudinal subintervals.'],
],[206,298])
p('The mixture is used for population evidence and model comparison. Its posterior is saved as a diagnostic; the primary classifier does <b>not</b> impose an additional 90% posterior or Gaussian-core cut. It is not the downstream sieve-plane GMM cleanup.')
p('Synthetic missing-band controls check that an empty geometric position is not automatically claimed as occupied. A raw replay background may still defeat discovery or population evidence; that behavior is an intended subject of the ifarm test.')

section('11. The training-count rule, defined precisely')
p('For each geometric band, take the training events that the geometry places near its center. Divide delta into intervals 1 percentage point wide by default. In each interval containing at least the configured minimum count, store the 0.5% and 99.5% quantiles of r. These are numerical bounds on where that band\'s training data were observed.')
p('For an event at continuous delta, interpolate the lower and upper r bounds between neighboring retained interval centers. Add padding of 0.04 in standardized r. The event must lie between these bounds. A gap may be bridged only when the neighboring centers are at most 2 percentage points apart; the half-bin edge convention is also applied. At least two retained intervals are required for a band envelope.')
table(['Setting','What changing it changes'],[
    ['coverage_min_events: 10 (primary), 5 (sensitivity)','Whether a 1%-delta interval contributes observed r bounds. This changes acceptance, not coefficients.'],
    ['coverage_bin_width: 1','The delta resolution and counting population of the training-envelope table.'],
    ['coverage_quantiles: [0.005, 0.995]','How far into the observed r tails the table reaches.'],
    ['coverage_padding: 0.04','Additional allowed r distance outside the interpolated quantile envelope.'],
    ['maximum_bridge: 2','Largest distance in delta between interval centers allowed to bridge a gap.'],
],[209,295])
p('An event on a visible density ridge can fail this rule when the random training subset has too few near-center events in its local delta interval, or when the event lies outside a tail quantile. The 3D cloud spans many delta intervals; nearby-looking points need not be in the same 1% interval.')
p('Changing only the count at prediction time leaves F, q, r, u, the proposed ridge and the population mixture fixed. New accepted events can reveal existing A/B assignment differences. Lowering the count does not create those fit differences. Bounds can also change non-monotonically as newly retained intervals alter interpolation.')

section('12. Rejection reasons and the X/Y join')
table(['Saved reason','Exact interpretation','First diagnostic'],[
    ['nonfinite','At least one FP coordinate or delta is not finite.','Replay branch values.'],
    ['delta_range','Delta is outside the configured half-open modeling range.','Metadata edges and arm acceptance.'],
    ['inverse_range','q is outside the outer half-index surfaces, or a locally increasing index direction cannot be inverted.','Measured FP plot with outer ridges; ordering checks.'],
    ['occupancy','The reporting interval has no eligible detected position for assignment.','Per-position population gains and counts.'],
    ['coverage','A nearest eligible position exists, but its event r/delta fails the training-envelope rule.','Near-center counts and r quantile bounds versus delta.'],
    ['accepted','All geometric, population and envelope conditions pass.','Physical naming and measured density.'],
],[90,237,177])
p('X and Y acceptance are separate. An X candidate can exist when Y is rejected, and conversely. The conventional downstream workflow can combine their masks using the preserved original replay entry. A complete hole identity requires both physical indices.')
p('The primary rejection reason records the first applicable limiting category according to the implementation\'s precedence. Other saved flags allow inspection when several conditions fail. A proposed nearest index, when one exists, is retained even for an envelope rejection.')
sub('Why write candidates directly')
p('The surface is evaluated at every original event delta, then conventional candidate fields are exported. Candidate tree names, branch types and summary tables remain compatible with the existing GMM consumers. Original replay entries are preserved; event order in a plotted sample is not an event identifier.')
p('These candidates are populations to inspect and possibly clean. Acceptance by a surface does not prove electron purity, a foil vertex, or absence of sieve scattering.')

section('13. What the seven-sample pilot found')
p('All seven samples were HMS at nominal momentum magnitude 6.667 GeV/c, near 12.49 degrees. They span delta from -10% to +10% within that setting; they do not test another nominal momentum. The total candidate pool was 421,666 events, with 84,338 common evaluation events.')
table(['Runs / foil z','Compared final pairs','Different pair','No new pair'],[
    ['1544 / 0 cm','41,604','211','17'],
    ['1540 / -8 cm','8,476','88','20'],
    ['1540 / +8 cm','7,270','5','25'],
    ['6985 / -3 cm','4,579','12','54'],
    ['6985 / +3 cm','4,897','19','61'],
    ['6983/6984 / -8 cm','4,483','70','53'],
    ['6983/6984 / +8 cm','3,524','11','85'],
    ['Total','74,833','416','315'],
],[155,120,111,118])
p('Fit A, minimum count 10: 74,417 events have the same new and old final pair; 416 differ; 315 final-labeled events lack a new pair. There are 8,701 newly assigned candidates absent from the old final selection. Of the 416 differences, 415 are X-only and one is Y-only. Overall compared-pair agreement is 99.44%. These are differences from an imperfect reference, not measured error rates.')
p('Independent A/B joined assignments differ on 144 of 83,311 events assigned by both. Differences concentrate at low delta. Agreement between two fits is reproducibility evidence; it does not by itself prove the physical identity.')
sub('The input restriction matters')
p('The central sample was a union of prior X/Y candidate selections. The other six samples were stricter intersections of prior X/Y candidates. The export applied no additional cleanup, but the old focal-plane cuts already shaped admission. The new raw-replay path removes that old-band admission and is therefore a substantive new test.')

section('14. Physical-coordinate comparisons')
fig('sieve_plane',height=270,caption='Run 1544 central foil. All common evaluation events in the two intervals, at their original reconstructed sieve positions. Color is the new X row; gray crosses lack a new X/Y pair. Nominal grid lines are guides. These positions were not moved, refocused or fitted by the surface algorithm.')
p('Elongation and overlap in a reconstructed sieve plane can coexist with distinct tracks in the two focal-plane spaces. The surface method assigns the row and column from those tracks; it does not optimize the optics matrix or make the sieve images narrower.')
sub('The difficult local comparison')
p('For runs 6983/6984, foil -8 cm, delta [-10,-8)%, 59 of 411 final-labeled events with a new X assignment differ from the old X row. All changes in this comparison are toward the next row. Many form a geometrically coherent sequence. That deserves physical inspection rather than tuning the method to reproduce every old label.')
p('The original reference was very good but was not perfect. Some historical components were deliberately split and later assigned together. The comparison uses final post-GMM event pairs, not intermediate cluster component IDs. Differences alone are not a failure criterion.')

section('14a. Inspecting a low-delta identity change')
fig('worst_low_delta_x_3d',height=330,caption='Identical 411 colored events in both 3D views, runs 6983/6984, foil -8 cm, -10% <= delta < -8%. Left: final post-GMM X rows. Right: new frozen surface X rows. Black outlines identify the 59 differing assignments; the faint background contains all 483 low-delta evaluation events. There is no Y-selection requirement in this comparison. Every event remains at the same measured FP coordinates and delta.')
p('Following a continuous measured ridge is more informative than counting agreement in the reconstructed sieve plane alone. This view makes the coherent changes toward the next row visible without asserting which assignment is physically correct.')
p('The new label is determined by the nearest eligible center of the global surface at each event\'s r and delta. The low-delta section is constrained by the same polynomial coefficients used elsewhere in delta. It is not a separate low-delta cluster fit, and it is not an event-by-event correction using reconstructed sieve proximity.')
p('The whole-family physical naming map is fixed across delta. Consequently a systematic change confined to the lowest delta can indicate a geometric interpretation difference, a boundary difference, or a problem in the historical local labeling. It cannot be repaired by quietly changing the global name only in that slice.')

section('15. What is ready, and what is still a test')
p('The farm branch contains a portable numerical engine, original-replay exporter, frozen-model and naming stages, conventional pre-GMM output adapter, reusable 3D inspection tools, and an explicit run guide. It does not import the local research output directories.')
table(['Established locally','Not established by the local pilot'],[
    ['Ordered geometry and reproducibility on seven preselected HMS samples.','Performance before the old X/Y band selections.'],
    ['Targeted numerical controls including absent populations, ordering, count-only acceptance changes and SHMS delta through +22%.','SHMS physics performance, another nominal momentum, and arbitrary initial reconstruction errors.'],
    ['Campaign format and original-entry compatibility are checked with synthetic ROOT inputs; the accompanying run guide records the checks actually completed.','Full-replay electron purity, sparse-tail provenance, foil-versus-sieve scattering, and improved focusing after a new SVD matrix.'],
],[252,252])
p('The HMS/SHMS folder prefix selects branch names, Cherenkov detector, acceptance and centered sieve geometry. It makes the software path compatible with both arms; it does not transfer measured HMS coefficients to SHMS. Each arm/foil learns its own geometry from its replay events.')
sub('A bounded first farm test')
p('Start with one familiar HMS rungroup. Inspect the stage-01 foil cuts and the new full pre-band cloud, then run fit, label and candidates. Compare physical names and rejection patterns with the prior campaign, especially the lowest delta. Repeat on a second foil or a SHMS campaign after this path is understood. No additional SVD run is required to inspect these candidate populations.')

section('16. Running on ifarm')
p('The exact environment checks, paths and stage outputs are documented in <b>JOINT_SURFACE_IFARM.md</b> at the repository root. Run from the checkout root in the normal HCANA/ROOT environment. Python needs NumPy and SciPy; ROOT I/O needs PyROOT and NumPy; 3D inspection needs Matplotlib with a GUI backend or Plotly for offline HTML.')
code('git fetch origin\ngit switch -c codex/joint-surface-ifarm \\\n  --track origin/codex/joint-surface-ifarm\n\npython3 setup_joint_surface_campaign.py \\\n  /path/to/stable/HMS_6p667GeV HMS_6p667GeV_joint_trial\nC=HMS_6p667GeV_joint_trial\npython3 run_joint_surface_campaign.py "$C" --plan\npython3 run_joint_surface_campaign.py "$C" ridge\npython3 run_joint_surface_campaign.py "$C" export\npython3 run_joint_surface_campaign.py "$C" fit\npython3 run_joint_surface_campaign.py "$C" label\npython3 run_joint_surface_campaign.py "$C" candidates')
p('If the branch already exists locally, switch to it and pull with --ff-only. Use a separate checkout if you want to preserve another active branch on ifarm. Source and destination campaigns must retain the same HMS_ or SHMS_ prefix. Setup copies input configuration and seed matrix, not numbered outputs or active old component overrides.')
sub('Inspect with the native viewer first')
code('python3 -m joint_surface.viewer "$C" --desktop')
p('Run in a graphical session, or use X forwarding from a machine with an X server. A headless batch process cannot display a rotatable desktop window. The HTML mode bundles its plotting library and can be opened locally in Firefox; see the run guide for the final command and troubleshooting.')

section('17. Tuning geometry when the result is poor')
table(['Symptom / setting','What to change and inspect'],[
    ['Bands tilted in the fixed view: theta_deg','Try a reviewed standardized angle, then refit in a fresh campaign. Do not choose separate angles per delta interval.'],
    ['Weak or doubled periodic proposal: min_q_frequency, max_q_frequency, phase_starts, phase_iterations','Inspect raw planes and q/r projections. These control the initialization search, not the hardware row count. Increasing the search does not establish that a weak mode is physical.'],
    ['Initializer not reproducible: minimum_validation_phase_power','Default 0.03. Prefer examining input/background and phase candidates before lowering the reproducibility gate.'],
    ['Excess curvature or delta dependence: families, high_order_penalty, extra_penalty','Restrict candidate families or strengthen penalties. Defaults 0.00002 for higher-order and 0.002 for extra delta terms. Keep the model-choice score and condition diagnostics.'],
    ['Model choice unstable: model_choice_tolerance, seed','Default tolerance 0.01 nat and fixed seed. Compare a separately seeded fit as an audit; do not select the seed by old-label agreement.'],
    ['Bad ridge centers: seed_radius, refine_radius, iterations, irls','Defaults 0.24/0.30 in index spacing, 4 outer iterations and 6 robust updates. These affect which near-center events fit geometry. More iterations cannot repair an incorrect periodicity automatically.'],
],[202,302])
p('Store shared settings under fit.common and axis-specific settings under fit.x or fit.y. Axis settings override common settings. Unknown engine keys are rejected. The exact defaults and selected settings are saved with each model.')
p('Any change affecting the frame, initializer, coefficients or population evidence requires refitting and redoing physical naming. Make a fresh campaign so outputs and candidate masks remain traceable.')

section('18. Tuning acceptance and physical naming')
table(['Problem','Relevant control','What remains fixed'],[
    ['Visible ridge events rejected by local training counts','coverage_min_events: try 5 versus primary 10, then inspect the recovered events by delta.','Coefficients and population fits can remain frozen.'],
    ['Events outside r-tail bounds','coverage_quantiles, coverage_padding','The reason is an observed training-envelope restriction, not proof that the ridge is absent.'],
    ['A gap in the training table','coverage_bin_width, maximum_bridge','Do not infer reliability across a large empty r/delta region merely by extending the bridge.'],
    ['An expected slot is not detected','occupancy_likelihood_gain, occupancy_min_effective_events, occupancy_min_coherent_intervals, occupancy_min_events','The surface predicts a location; measured density must still support a population.'],
    ['All physical rows/columns shifted or reversed','global_maps orientation and origin after geometry review','No center, boundary or latent assignment changes.'],
    ['A single band is merged or doubled','Initialization and family fit, not event-specific naming patches','A global map cannot repair a geometric count/index error.'],
],[136,184,184])
p('The first branch release deliberately protects written stages from replacement. To change settings after a stage has been run, make a fresh campaign and repeat the necessary stages. The run guide explains a count-only sensitivity rerun and what it costs. No knob is a target for maximizing historical-label agreement.')
p('The defaults do not include a separate tail cleanup based on electron/vertex provenance. That study remains open. The density cleanup in another branch is also separate. These should be explicit later selections with saved masks, not hidden inside the surface fit.')

section('19. Mathematical taxonomy and reading onramp')
table(['Technique family','Its precise role here'],[
    ['Unsupervised latent-index estimation','Periodic density discovers temporary integer positions without physical training labels. This is structured geometric discovery.'],
    ['Conditional parametric regression','F predicts q from index, r and delta with a finite polynomial basis. Nonlinear coordinates do not make it nonlinear in coefficients.'],
    ['Alternating estimation','Geometry-inferred indices and coefficients are updated in turn. This resembles latent-variable methods, but is not a proof of maximum-likelihood EM convergence.'],
    ['Robust M-estimation / iteratively reweighted least squares','Reduce the influence of large q residuals in center fits; the smooth weighting is related to robust losses.'],
    ['Generalized ridge regularization','A diagonal penalty shrinks different coefficient groups by different amounts. It limits variance and ill-conditioning.'],
    ['Shape-constrained regression','Positive change with ordered index prevents local order reversal over the checked domain.'],
    ['Constrained mixture density estimation','Training-only population tests and normalized conditional density scores; this is distinct from downstream sieve-plane GMM cleanup.'],
    ['Geometric classification with rejection','Nearest eligible sheet plus explicit inversion, population and training-envelope requirements.'],
],[190,314])
p('Principal curves and surfaces offer a useful geometric analogy, but this is not the self-consistent principal-curve algorithm. It is not a general tensor-product spline, neural network, or joint X-Y density model. The "ridge" in ridge regression denotes a coefficient regularization method, not a detector density ridge.')

section('20. Canonical sources and next decisions')
p('Read the regression and model-selection material first, then the optimization and latent-variable references. These sources explain method families; none is a publication validating this specific Hall C implementation.','small')
for text in [
    '<b>Hastie, Tibshirani and Friedman (2009).</b> <i>The Elements of Statistical Learning</i>, 2nd ed., Springer. Chapters 3, 5, 7 and 14: basis regression, regularization, independent model assessment and unsupervised learning. <link href="https://hastie.su.domains/ElemStatLearn/main.html">Author book page</link>.',
    '<b>Hoerl and Kennard (1970).</b> "Ridge Regression: Biased Estimation for Nonorthogonal Problems," <i>Technometrics</i> 12(1), 55-67. DOI: 10.1080/00401706.1970.10488634. Generalized diagonal shrinkage is the relevant connection; it cannot create missing covariate information.',
    '<b>Huber (1964).</b> "Robust Estimation of a Location Parameter," <i>Annals of Mathematical Statistics</i> 35(1), 73-101. DOI: 10.1214/aoms/1177703732. A foundation for reduced influence of large residuals; the implemented smooth weights are not the literal piecewise Huber loss.',
    '<b>Boyd and Vandenberghe (2004).</b> <i>Convex Optimization</i>, Cambridge. Sections 4.4 and 6.3-6.5: constrained quadratic fitting, regularization and robust approximation. <link href="https://web.stanford.edu/~boyd/cvxbook/">Official freely available book</link>. The coefficient subproblem does not make global periodic discovery convex.',
    '<b>Dempster, Laird and Rubin (1977).</b> "Maximum Likelihood from Incomplete Data via the EM Algorithm," <i>JRSS B</i> 39(1), 1-22. DOI: 10.1111/j.2517-6161.1977.tb01600.x. Read for the mixture-fitting connection; alternating hard index updates here are not this theorem.',
    '<b>Hastie and Stuetzle (1989).</b> "Principal Curves," <i>JASA</i> 84(406), 502-516. DOI: 10.1080/01621459.1989.10478797. A useful geometric comparison, not the algorithm implemented here.',
    '<b>Silverman (1986).</b> <i>Density Estimation for Statistics and Data Analysis</i>, Chapman and Hall. Read for the effects of smoothing and projection on modes. A projected peak is not automatically a physical sieve identity.',
]:p(text,'small')
sub('What should proceed to SVD')
p('A production fit should use one explicitly chosen set of physical X/Y event masks, followed by the chosen cleanup and the existing fit-ntuple construction. Keep ambiguous or rejected events identifiable. Reproducibility audits inform that choice; they do not supply duplicate conflicting rows to the matrix solver. This branch stops at pre-GMM candidates so that the farm test and any later cleanup remain reviewable.')
p('Next: inspect full replay on ifarm, confirm global physical naming, and compare raw-plane density and rejection patterns. Then return to the SVD optimization with the accepted candidate/cleanup policy stated explicitly.','small')


def footer(canvas,doc):
    canvas.setStrokeColor(colors.HexColor('#cad4dc'))
    canvas.line(54,43,558,43)
    canvas.setFont('Helvetica',8)
    canvas.setFillColor(MUTED)
    canvas.drawString(54,30,'Continuous band surfaces | Farm-test method guide | 27 October 2026')
    canvas.drawRightString(558,30,str(doc.page))


if __name__=='__main__':
    document=SimpleDocTemplate(str(OUT),pagesize=letter,leftMargin=54,rightMargin=54,
                              topMargin=43,bottomMargin=55,title='Continuous band surfaces: HMS and SHMS method guide',
                              author='Hall C Optics provisional study')
    document.build(story,onFirstPage=footer,onLaterPages=footer)
    print(OUT)
