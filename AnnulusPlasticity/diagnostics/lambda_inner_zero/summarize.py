"""Summarize direct constitutive-branch traces from instrumented SBM solves."""
import csv
import json
from pathlib import Path
import re

OUT = Path(__file__).resolve().parent
FIELDS = ('time', 'iteration', 'condition', 'D_s_elastic_distance',
          'D_t_elastic_distance', 'alpha_s_old', 'alpha_t_old',
          'alpha_s_trial', 'alpha_t_trial', 'x_s', 'y_s', 'x_t', 'y_t')
D_TOL = 1e-8
ALPHA_TOL = 1e-12


def main():
    summaries = []
    all_iterations = []
    for path in sorted(OUT.glob('lambda_*/*/n*/result.json')):
        result = json.loads(path.read_text())
        groups = {}
        with (path.parent/'run.log').open() as stream:
            for line in stream:
                match = re.search(r'REGIME_PAIR\s+(.+)', line)
                if not match:
                    continue
                values = match.group(1).split()
                if len(values) != len(FIELDS):
                    raise ValueError(line)
                row = {name: float(value) for name, value in zip(FIELDS, values)}
                for name in ('iteration', 'condition'):
                    row[name] = int(row[name])
                row['active_s'] = row['D_s_elastic_distance'] > D_TOL
                row['active_t'] = row['D_t_elastic_distance'] > D_TOL
                row['plasticized_s'] = row['alpha_s_trial'] > ALPHA_TOL
                row['plasticized_t'] = row['alpha_t_trial'] > ALPHA_TOL
                row['committed_plasticized_s'] = row['alpha_s_old'] > ALPHA_TOL
                row['committed_plasticized_t'] = row['alpha_t_old'] > ALPHA_TOL
                key = (row['time'], row['iteration'])
                groups.setdefault(key, {})[row['condition']] = row
        if not groups:
            raise ValueError(f'No constitutive traces in {path.parent}')
        records = []
        for (time, iteration), conditions in sorted(groups.items(), key=lambda x: (x[0][0], x[0][1] < 0, abs(x[0][1]))):
            rows = list(conditions.values())
            record = dict(lambda_inner=result['lambda_inner'], spans=result['spans'], strategy=result['strategy'],
                time=time, iteration=abs(iteration), final_snapshot=iteration < 0,
                points=len(rows),
                elastic_s_plastic_t=sum(not r['active_s'] and r['active_t'] for r in rows),
                plastic_s_elastic_t=sum(r['active_s'] and not r['active_t'] for r in rows),
                both_elastic=sum(not r['active_s'] and not r['active_t'] for r in rows),
                both_active_plastic=sum(r['active_s'] and r['active_t'] for r in rows),
                accumulated_mismatch=sum(r['plasticized_s'] != r['plasticized_t'] for r in rows),
                committed_mismatch=sum(r['committed_plasticized_s'] != r['committed_plasticized_t'] for r in rows),
                active_mismatch_both_with_plastic_history=sum(
                    r['active_s'] != r['active_t'] and r['plasticized_s'] and r['plasticized_t'] for r in rows))
            record['active_mismatch'] = record['elastic_s_plastic_t'] + record['plastic_s_elastic_t']
            records.append(record)
        terminal_time = result['time']
        terminal = next(r for r in records if r['final_snapshot'] and r['time'] == terminal_time)
        summary = dict(lambda_inner=result['lambda_inner'], geometry=result['geometry'], L2_relative=result.get('L2_relative'), H1_bulk_relative=result.get('H1_bulk_relative'), spans=result['spans'], strategy=result['strategy'],
            status=result['status'], time=terminal_time,
            max_active_mismatch=max(r['active_mismatch'] for r in records),
            max_accumulated_mismatch=max(r['accumulated_mismatch'] for r in records),
            terminal=terminal,
            final_steps=[r for r in records if r['final_snapshot']])
        summaries.append(summary)
        all_iterations.extend(records)
        print(json.dumps({k:v for k,v in summary.items() if k != 'final_steps'}), flush=True)
        with (path.parent/'pairs.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(next(iter(next(iter(groups.values())).values()))))
            writer.writeheader()
            for conditions in groups.values():
                writer.writerows(conditions.values())
        with (path.parent/'iterations.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    (OUT/'summary.json').write_text(json.dumps(summaries, indent=2)+'\n')
    with (OUT/'iterations.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_iterations[0]))
        writer.writeheader()
        writer.writerows(all_iterations)


if __name__ == '__main__':
    main()
