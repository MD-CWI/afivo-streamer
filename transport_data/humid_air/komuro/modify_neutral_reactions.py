#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv('neutral_reactions.txt', dtype=str)

print('number,reaction,rate,comment')

for index, row in df.iterrows():
    row['reaction'] = row['reaction'].strip()
    row['rate'] = row['rate'].strip()
    if float(row['beta']) != 0:
        row['rate'] += f'*(Tg/300)**{row["beta"].strip()}'
    if float(row['Ea']) != 0:
        row['rate'] += f'*exp(-{row["Ea"].strip()}/Tg)'
        row['rate'] = row['rate'].replace('--', '')
    print(f"{row['number']},{row['reaction']},{row['rate']},{row['comment']}")