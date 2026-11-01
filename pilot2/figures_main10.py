"""Readable main-text figures for the ten-page revision (31 October 2026).

python -m pilot2.figures_main10 --out paper/figures

Reuse the published figure_data.csv, interpolation, seed aggregation and saturation
reader from figures_privacy. Only layout and labeling change. The single-copy
borrowed-host control is shown only in its matching left panel, not reused as a
sixteen-copy donor-tier control. No experiment results are rewritten.
"""
import argparse
from pathlib import Path
from . import figures_privacy as f
import matplotlib.pyplot as plt


def save(fig, out, name):
    # Exact 6.5-inch media box: 10.2 pt fonts stay >=10 pt at manuscript width.
    fig.savefig(out / (name + '.pdf'), bbox_inches=None, pad_inches=0)
    fig.savefig(out / (name + '.png'), bbox_inches=None, pad_inches=0, dpi=220)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=Path('paper/figures'))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = f.load(Path('results/figure_data.csv'))
    plt.rcParams.update({'font.size': 10.2, 'axes.titlesize': 10.2,
                         'axes.labelsize': 10.2, 'legend.fontsize': 10.2,
                         'xtick.labelsize': 10.2, 'ytick.labelsize': 10.2})

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.8), sharey=True)
    fig.subplots_adjust(left=.09, right=.985, top=.90, bottom=.39, wspace=.16)
    runs = {f'dnabert2bpe_human_s{s}' for s in (0, 1, 2)}
    for ax, tier, title in zip(axes, ('1', '16'), ('Single copy', 'Sixteen copies')):
        for host, field, color, label, marker in (
            ('train', 'rank1_rate', f.C_ID, 'Membership: own', 'o'),
            ('fresh', 'rank1_rate', f.C_ID2, 'Membership: fresh', 'd'),
            ('train', 'exact_rate', f.C_VERB, 'Prefix recovery: own', 's')):
            f._plot_band(ax, f._band(rows, 'bpe_human', runs, tier, host, field),
                         color, label, marker)
        if tier == '1':
            for host, label, marker in (
                ('train', 'Control: borrowed', 'v'),
                ('fresh', 'Control: fresh', '^')):
                f._plot_band(ax, f._band(rows, 'bpe_human', runs, '0', host, 'rank1_rate'),
                             f.C_REF, label, marker)
        ax.axhline(1 / 101, color=f.C_REF, lw=.7, ls=':')
        ax.set_title(title)
        ax.set_xlim(0, 160); ax.set_ylim(-.045, 1.04)
        ax.set_xticks([15, 50, 100, 150], ['15', '50', '100', 'final'])
        ax.set_xlabel('Fine-tuning round', labelpad=2)
        f._style(ax)
    axes[0].set_ylabel('Rate'); axes[0].set_yticks([0, .5, 1], ['0', '0.5', '1.0'])
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc='lower center', bbox_to_anchor=(.51, .005),
               ncol=2, frameon=False, handlelength=1.7, columnspacing=1.4,
               labelspacing=.25, borderaxespad=0)
    save(fig, args.out, 'F6_ident')

    main_data, repl = f.matched(rows, 'kcurve300'), f.matched(rows, 'kcurve')
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.9))
    fig.subplots_adjust(left=.09, right=.985, top=.91, bottom=.38, wspace=.40)
    ax, bx = axes
    toks = [t for t in f.TOKS if t in main_data]
    x = list(range(len(toks)))
    colors = ['#9CAAB9', '#768EAB', '#4F6B90', '#1E4F8A']
    for ti, target in enumerate(f.TARGETS):
        valid = [(i, main_data[t][0][ti]) for i,t in enumerate(toks)
                 if main_data[t][0][ti] is not None]
        ax.plot([i for i,y in valid], [y for i,y in valid], '-o',
                ms=3.5, lw=1.1, color=colors[ti], label=f'{target:.2f} bits/nt')
    valid = [(i,repl[t][0][-1]) for i,t in enumerate(toks)
             if t in repl and repl[t][0][-1] is not None]
    ax.plot([i for i,y in valid], [y for i,y in valid], 'o', mfc='none',
            mec=f.INK, ms=6, mew=.8, label='288 nt, at 2.40')
    ax.set_title('Matched-loss recovery')
    ax.set_ylabel('Exact recovery')
    ax.set_xticks(x, [f'{i+1}\n{f.SLOTS[t]:g}' for i,t in enumerate(toks)])
    ax.set_xlabel('Token width k / nominal slots', labelpad=2)
    ax.set_ylim(-.04, 1.04); ax.set_yticks([0, .5, 1], ['0', '50%', '100%'])
    ax.annotate('Unreached\nlosses', (0,.1), color=f.C_REF, fontsize=10.2)
    f._style(ax)
    for table, marker, color, fill, label in (
        (main_data, 's', f.C_ID, f.C_ID, '300 nt'),
        (repl, 'o', f.INK, 'none', '288 nt')):
        valid = [(i, table[t][1]) for i,t in enumerate(toks)
                 if t in table and table[t][1] is not None]
        bx.plot([i for i,y in valid], [y for i,y in valid], '-'+marker,
                color=color, mfc=fill, ms=4, lw=1, label=label)
    bx.set_title('Full-recovery loss')
    bx.set_ylabel('Held-out bits/nt')
    bx.set_xticks(x, [str(i+1) for i in x]); bx.set_xlabel('Token width k', labelpad=2)
    bx.set_xlim(-.2,5.3); bx.set_ylim(2.3,4.5)
    bx.set_yticks([2.5,3.5,4.5]); bx.text(0,2.45,'Not reached',color=f.C_REF,fontsize=10.2)
    bx.legend(frameon=False,loc='upper right',handlelength=1.3)
    f._style(bx)
    h,lab=ax.get_legend_handles_labels()
    fig.legend(h,lab,loc='lower center',bbox_to_anchor=(.5,.005),ncol=3,
               frameon=False,handlelength=1.4,columnspacing=1.2,labelspacing=.3)
    save(fig,args.out,'F5_slots')


if __name__ == '__main__':
    main()
