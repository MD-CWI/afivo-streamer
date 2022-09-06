We show a reaction list in Komuro_reactions.csv for N2-O2 mixtures from Ono and Komuro 2020 (https://doi.org/10.1088/1361-6463/ab4e65).

We construct a new reaction list in Baohong_reactions.csv for N2-O2 mixtures, which are primarily from Kossyi 1992 (https://doi.org/10.1088/0963-0252/1/3/011),  Fresnet 2002 (https://doi.org/10.1088/0963-0252/11/2/305), Gordillo 2008 (https://doi.org/10.1088/0022-3727/41/23/234016), Popov 2011 (https://doi.org/10.1088/0022-3727/44/28/285201), Panchesnyi 2013 (https://doi.org/10.1088/0022-3727/46/15/155201) and Ono and Komuro 2020 (https://doi.org/10.1088/1361-6463/ab4e65).

The csv format of reaction list can be converted to the format of reactions used in Afivo-streamer via "chemistry_reaction_parser.py":

python chemistry_reaction_parser.py Komuro_reactions.csv -convert_tex Komuro_converted.csv 

python chemistry_reaction_parser.py Komuro_converted.csv
