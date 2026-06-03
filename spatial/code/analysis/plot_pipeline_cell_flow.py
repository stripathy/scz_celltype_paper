"""
Pipeline cell flow river/Sankey plot.

Traces cells through each QC stage from raw to final hybrid_qc_pass,
showing where cells are flagged, rescued, and excluded.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Load data ────────────────────────────────────────────────────────
df = pd.read_csv('output/pipeline_cell_flow.csv')

# ── Compute pipeline-wide totals ─────────────────────────────────────
T = {}
T['total'] = df['n_total'].sum()
T['qc_pass'] = df['qc_pass'].sum()
T['qc_fail_other'] = df['qc_fail'].sum() - df['fail_total_high'].sum()
T['fail_total_high'] = df['fail_total_high'].sum()
T['high_umi_only'] = df['high_umi_only'].sum()
T['high_umi_multi_fail'] = T['fail_total_high'] - T['high_umi_only']
T['low_margin'] = df['low_margin_fail'].sum()
T['doublet_suspect'] = df['doublet_suspect'].sum()
T['glut_gaba'] = df['glut_gaba_doublet'].sum()
T['gaba_gaba'] = df['gaba_gaba_doublet'].sum()
T['resolved'] = df['nuc_resolved'].sum()
T['persistent'] = df['nuc_persistent'].sum()
T['nuclear_only'] = df['nuc_nuclear_only'].sum()
T['insufficient'] = df['nuc_insufficient'].sum()
T['hybrid_qc_pass'] = df['hybrid_qc_pass'].sum()

# ── FIGURE 1: Three-panel overview ───────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(16, 22),
                         gridspec_kw={'height_ratios': [2.5, 1.5, 2]})

# ============================================================
# PANEL A: Delta bar chart (gains and losses at each step)
# ============================================================
ax = axes[0]

# Each event is a delta from the pipeline
events = [
    ('QC fail:\ncontrol probes\n& low genes', -T['qc_fail_other'], '#CC6677'),
    ('QC fail:\nhigh UMI\n(multi-fail)', -T['high_umi_multi_fail'], '#CC6677'),
    ('QC fail:\nhigh UMI only\n(→ rescued later)', -T['high_umi_only'], '#DDAA33'),
    ('Low corr\nmargin\n(bottom 5%)', -T['low_margin'], '#CC6677'),
    ('Doublet:\nGlut+GABA', -T['glut_gaba'], '#CC6677'),
    ('Doublet:\nGABA+GABA', -T['gaba_gaba'], '#CC6677'),
    ('High-UMI\nrescued\n(Step 04)', +T['high_umi_only'], '#228833'),
    ('Doublets\nresolved\n(nuclear clear)', +T['resolved'], '#228833'),
    ('Persistent\ndoublets', -T['persistent'], '#CC6677'),
    ('Nuclear-only\ndoublets', -T['nuclear_only'], '#CC6677'),
    ('Insufficient\nevidence', -T['insufficient'], '#882255'),
]

x_pos = np.arange(len(events))
vals = [e[1] for e in events]
colors = [e[2] for e in events]
labels = [e[0] for e in events]

bars = ax.bar(x_pos, vals, color=colors, width=0.65, edgecolor='white', linewidth=0.5)

# Add value labels
for i, (xp, v) in enumerate(zip(x_pos, vals)):
    if abs(v) < 500:
        # Too small to label inside bar
        yoff = -1500 if v < 0 else 1500
        ax.text(xp, v + yoff, f'{v:+,}', ha='center', va='center',
                fontsize=11, fontweight='bold', color=colors[i])
    else:
        ax.text(xp, v/2, f'{v:+,}\n({abs(v)/T["total"]*100:.2f}%)',
                ha='center', va='center', fontsize=11, fontweight='bold',
                color='white' if abs(v) > 3000 else colors[i])

ax.axhline(0, color='black', linewidth=0.8)
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel('Cell count change', fontsize=14)
ax.set_title(f'Pipeline Cell Flow: Gains and Losses at Each QC Step\n'
             f'(Start: {T["total"]:,} → Final: {T["hybrid_qc_pass"]:,} = '
             f'{T["hybrid_qc_pass"]/T["total"]*100:.1f}% retained, '
             f'{T["total"]-T["hybrid_qc_pass"]:,} lost)',
             fontsize=16, fontweight='bold')
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e3:+.0f}K' if abs(x) >= 1000 else f'{x:+.0f}'))

# Add step annotations
ax.axvspan(-0.5, 2.5, alpha=0.06, color='blue', zorder=0)
ax.axvspan(2.5, 5.5, alpha=0.06, color='orange', zorder=0)
ax.axvspan(5.5, 10.5, alpha=0.06, color='green', zorder=0)
ax.text(1, ax.get_ylim()[0]*0.92, 'Step 01: Spatial QC', ha='center',
        fontsize=12, fontstyle='italic', color='#333')
ax.text(4, ax.get_ylim()[0]*0.92, 'Step 02b: Corr Classifier', ha='center',
        fontsize=12, fontstyle='italic', color='#333')
ax.text(8, ax.get_ylim()[0]*0.92, 'Step 04: Nuclear Resolution', ha='center',
        fontsize=12, fontstyle='italic', color='#333')

legend_elements = [
    mpatches.Patch(facecolor='#CC6677', label='Cells excluded'),
    mpatches.Patch(facecolor='#DDAA33', label='Pending rescue (high-UMI only)'),
    mpatches.Patch(facecolor='#228833', label='Cells rescued'),
    mpatches.Patch(facecolor='#882255', label='Insufficient evidence'),
]
ax.legend(handles=legend_elements, loc='lower left', fontsize=12)

# ============================================================
# PANEL B: Running total through pipeline
# ============================================================
ax1 = axes[1]

# Build running total
checkpoints = [
    ('Raw\ncells', T['total']),
    ('After\nStep 01\n(qc_pass)', T['qc_pass']),
    ('After\nStep 02b\n(corr_qc_pass)', T['qc_pass'] - T['low_margin'] - T['doublet_suspect']),
    ('After Step 04\n(hybrid_qc_pass)', T['hybrid_qc_pass']),
]

cp_x = np.arange(len(checkpoints))
cp_vals = [c[1] for c in checkpoints]
cp_labels = [c[0] for c in checkpoints]

ax1.plot(cp_x, cp_vals, 'o-', color='#4477AA', linewidth=3, markersize=12, zorder=5)

for i, (xp, v) in enumerate(zip(cp_x, cp_vals)):
    pct = v / T['total'] * 100
    ax1.annotate(f'{v:,}\n({pct:.1f}%)', (xp, v),
                textcoords='offset points', xytext=(0, 18),
                ha='center', fontsize=13, fontweight='bold')

# Shade the rescue region
ax1.fill_between([1.5, 3.5],
                 [T['qc_pass'] - T['low_margin'] - T['doublet_suspect']]*2,
                 [T['hybrid_qc_pass']]*2,
                 alpha=0.15, color='#228833', zorder=0)
ax1.annotate(f'Net rescue:\n+{T["hybrid_qc_pass"] - (T["qc_pass"] - T["low_margin"] - T["doublet_suspect"]):,}',
             xy=(2.5, (T['hybrid_qc_pass'] + T['qc_pass'] - T['low_margin'] - T['doublet_suspect'])/2),
             fontsize=11, ha='center', va='center', color='#228833', fontweight='bold')

ax1.set_xticks(cp_x)
ax1.set_xticklabels(cp_labels, fontsize=12)
ax1.set_ylabel('Cells passing', fontsize=14)
ax1.set_title('Running Cell Count Through Pipeline', fontsize=16, fontweight='bold')
ax1.set_ylim(T['total'] * 0.92, T['total'] * 1.03)
ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.2f}M'))
ax1.grid(axis='y', alpha=0.3)

# ============================================================
# PANEL C: Per-sample breakdown
# ============================================================
ax2 = axes[2]

df_sorted = df.sort_values('n_total', ascending=True).reset_index(drop=True)
y_pos = np.arange(len(df_sorted))
bar_height = 0.7

# Final pass
ax2.barh(y_pos, df_sorted['hybrid_qc_pass'], height=bar_height,
         color='#4477AA', label='hybrid_qc_pass', alpha=0.9)

# QC fail (non-high-UMI, non-rescuable)
qc_fail_other = df_sorted['qc_fail'] - df_sorted['high_umi_only']
# Among qc_pass cells, those excluded by doublet/margin but not rescued
net_post_qc_excluded = df_sorted['n_total'] - df_sorted['hybrid_qc_pass'] - qc_fail_other

# Stack: high-UMI rescued (already in hybrid_qc_pass) shown separately
ax2.barh(y_pos, qc_fail_other.clip(lower=0), height=bar_height,
         left=df_sorted['hybrid_qc_pass'], color='#CC6677',
         label='QC fail (not rescued)', alpha=0.9)

ax2.barh(y_pos, net_post_qc_excluded.clip(lower=0), height=bar_height,
         left=df_sorted['hybrid_qc_pass'] + qc_fail_other.clip(lower=0),
         color='#DDAA33', label='Doublet/margin excluded (net)', alpha=0.9)

ax2.set_yticks(y_pos)
ax2.set_yticklabels(df_sorted['sample'], fontsize=11)
ax2.set_xlabel('Number of cells', fontsize=14)
ax2.set_title('Per-Sample Cell Fate (hybrid_qc_pass %)', fontsize=16, fontweight='bold')
ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e3:.0f}K'))
ax2.legend(loc='lower right', fontsize=11)

for idx, row in df_sorted.iterrows():
    pct = row['hybrid_qc_pass'] / row['n_total'] * 100
    ax2.text(row['n_total'] + 500, y_pos[idx], f'{pct:.1f}%',
             va='center', fontsize=10, color='#333333', fontweight='bold')

plt.tight_layout()
plt.savefig('output/presentation/pipeline_cell_flow.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: output/presentation/pipeline_cell_flow.png')

# ── FIGURE 2: Per-sample detail table ────────────────────────────────
fig2, ax3 = plt.subplots(figsize=(18, 10))
ax3.axis('off')

# Build table data
table_data = []
columns = ['Sample', 'Total', 'QC fail\n(probes/genes)', 'High-UMI\nonly fail',
           'qc_pass', 'Low\nmargin', 'Doublet\nsuspect', 'Resolved\n(rescued)',
           'Persistent', 'Nuclear\nonly', 'hybrid_qc\n_pass', 'Pass %']

for _, row in df.sort_values('sample').iterrows():
    qc_fail_non_high = row['qc_fail'] - row['fail_total_high']
    table_data.append([
        row['sample'],
        f"{int(row['n_total']):,}",
        f"{int(qc_fail_non_high):,}",
        f"{int(row['high_umi_only']):,}",
        f"{int(row['qc_pass']):,}",
        f"{int(row['low_margin_fail']):,}",
        f"{int(row['doublet_suspect']):,}",
        f"{int(row['nuc_resolved']):,}",
        f"{int(row['nuc_persistent']):,}",
        f"{int(row['nuc_nuclear_only']):,}",
        f"{int(row['hybrid_qc_pass']):,}",
        f"{row['hybrid_qc_pass']/row['n_total']*100:.1f}%"
    ])

# Add totals row
table_data.append([
    'TOTAL',
    f"{T['total']:,}",
    f"{T['qc_fail_other']:,}",
    f"{T['high_umi_only']:,}",
    f"{T['qc_pass']:,}",
    f"{T['low_margin']:,}",
    f"{T['doublet_suspect']:,}",
    f"{T['resolved']:,}",
    f"{T['persistent']:,}",
    f"{T['nuclear_only']:,}",
    f"{T['hybrid_qc_pass']:,}",
    f"{T['hybrid_qc_pass']/T['total']*100:.1f}%"
])

table = ax3.table(cellText=table_data, colLabels=columns,
                  loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.0, 1.4)

# Color headers
for j in range(len(columns)):
    table[0, j].set_facecolor('#4477AA')
    table[0, j].set_text_props(color='white', fontweight='bold')

# Color totals row
for j in range(len(columns)):
    table[len(table_data), j].set_facecolor('#E8E8E8')
    table[len(table_data), j].set_text_props(fontweight='bold')

# Color loss columns
loss_cols = [2, 5, 6, 8, 9]
rescue_cols = [3, 7]
for i in range(1, len(table_data) + 1):
    for j in loss_cols:
        table[i, j].set_facecolor('#FFF0F0')
    for j in rescue_cols:
        table[i, j].set_facecolor('#F0FFF0')

ax3.set_title('Per-Sample Pipeline Cell Flow Detail', fontsize=16, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig('output/presentation/pipeline_cell_flow_table.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: output/presentation/pipeline_cell_flow_table.png')

# ============================================================
# Print summary table
# ============================================================
print()
print('=' * 80)
print('PIPELINE CELL FLOW SUMMARY')
print('=' * 80)
print()
print(f'{"Stage":<40} {"Cells":>10} {"% of total":>10} {"Net Δ":>10}')
print('-' * 70)
print(f'{"Raw cells (Step 00)":<40} {T["total"]:>10,} {"100.0%":>10} {"":>10}')
print()
print(f'  Step 01 QC failures:')
print(f'    {"fail_neg_probe (>P99)":<36} {-5399:>10,} {5399/T["total"]*100:>9.2f}%')
print(f'    {"fail_neg_codeword (>P99)":<36} {-6176:>10,} {6176/T["total"]*100:>9.2f}%')
print(f'    {"fail_unassigned (>P99)":<36} {-7871:>10,} {7871/T["total"]*100:>9.2f}%')
print(f'    {"fail_n_genes_low (<med-5MAD)":<36} {-6:>10,} {6/T["total"]*100:>9.4f}%')
print(f'    {"fail_total_counts_low":<36} {0:>10,} {"0.00%":>10}')
print(f'    {"fail_total_counts_high (>med+5MAD)":<36} {-T["fail_total_high"]:>10,} {T["fail_total_high"]/T["total"]*100:>9.2f}%')
print(f'      {"of which high-UMI-only (rescuable)":<34} {T["high_umi_only"]:>10,} {T["high_umi_only"]/T["total"]*100:>9.2f}%')
print(f'      {"of which multi-fail (not rescuable)":<34} {T["high_umi_multi_fail"]:>10,} {T["high_umi_multi_fail"]/T["total"]*100:>9.2f}%')
print(f'  {"(Note: cells can fail multiple criteria)":<40}')
print()
print(f'{"After Step 01 (qc_pass)":<40} {T["qc_pass"]:>10,} {T["qc_pass"]/T["total"]*100:>9.2f}% {T["qc_pass"]-T["total"]:>+10,}')
print()
print(f'  Step 02b additional flags (among qc_pass):')
print(f'    {"Low margin (bottom 5% per sample)":<36} {-T["low_margin"]:>10,} {T["low_margin"]/T["total"]*100:>9.2f}%')
print(f'    {"Doublet: Glut+GABA":<36} {-T["glut_gaba"]:>10,} {T["glut_gaba"]/T["total"]*100:>9.2f}%')
print(f'    {"Doublet: GABA+GABA":<36} {-T["gaba_gaba"]:>10,} {T["gaba_gaba"]/T["total"]*100:>9.2f}%')
print()
corr_pass = T['qc_pass'] - T['low_margin'] - T['doublet_suspect']
print(f'{"After Step 02b (corr_qc_pass)":<40} {corr_pass:>10,} {corr_pass/T["total"]*100:>9.2f}% {corr_pass-T["qc_pass"]:>+10,}')
print()
print(f'  Step 04 nuclear doublet resolution:')
print(f'    {"High-UMI cells rescued":<36} {+T["high_umi_only"]:>+10,} {T["high_umi_only"]/T["total"]*100:>9.2f}%')
print(f'    {"Doublets resolved (nuclear clear)":<36} {+T["resolved"]:>+10,} {T["resolved"]/T["total"]*100:>9.2f}%')
print(f'    {"Persistent doublets (excluded)":<36} {-T["persistent"]:>10,} {T["persistent"]/T["total"]*100:>9.2f}%')
print(f'    {"Nuclear-only doublets (excluded)":<36} {-T["nuclear_only"]:>10,} {T["nuclear_only"]/T["total"]*100:>9.2f}%')
print(f'    {"Insufficient evidence (excluded)":<36} {-T["insufficient"]:>10,} {T["insufficient"]/T["total"]*100:>9.2f}%')
print()
print(f'{"FINAL (hybrid_qc_pass)":<40} {T["hybrid_qc_pass"]:>10,} {T["hybrid_qc_pass"]/T["total"]*100:>9.2f}% {T["hybrid_qc_pass"]-T["total"]:>+10,}')
print()
print(f'Net loss: {T["total"]-T["hybrid_qc_pass"]:,} cells ({(T["total"]-T["hybrid_qc_pass"])/T["total"]*100:.2f}%)')
print(f'Resolution rate: {T["resolved"]}/{T["resolved"]+T["persistent"]+T["insufficient"]+T["nuclear_only"]} doublet suspects '
      f'({T["resolved"]/(T["resolved"]+T["persistent"]+T["insufficient"]+T["nuclear_only"])*100:.1f}% rescued)')
