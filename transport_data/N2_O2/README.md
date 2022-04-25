TODO: briefly describe the chemistry dataset collected here
We collected the reaction for N2-O2 mixtures from Ono and Komuro 2020 (https://doi.org/10.1088/1361-6463/ab4e65), as listed in Komuro_reactions.csv
We created the reaction list for N2-O2 mixtures from Ono and Komuro 2020 (https://doi.org/10.1088/1361-6463/ab4e65), Panchesnyi 2013 (https://doi.org/10.1088/0022-3727/46/15/155201) and Kossyi 1992 (https://doi.org/10.1088/0963-0252/1/3/011), as listed in Baohong_reactions.csv
The latex format of the reaction list can be converted to the format of reactions used in Afivo-streamer via chemistry_reaction_parser.py:
python chemistry_reaction_parser.py Komuro_reactions.csv -convert_tex Komuro_converted.csv 
python chemistry_reaction_parser.py Komuro_converted.csv
