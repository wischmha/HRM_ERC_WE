#!/bin/sh
#
# Hans-Aloys Wischmann, 2024-07-21
# Ubuntu, or Ubuntu shell (app) on Windows (WSL 2)
# 
mkdir ../Downloads
mkdir ../Downloads/Population
mkdir ../Downloads/Geography
mkdir ../Downloads/Temperature
cd    ../Downloads
echo "Starting downloads:" `date` | tee Downloads.txt

cd Geography

wget https://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/JRC_GRID_2018.zip -O popd.zip
wget https://gisco-services.ec.europa.eu/distribution/v2/nuts/shp/NUTS_RG_01M_2021_4326.shp.zip -O nuts.zip
wget https://gisco-services.ec.europa.eu/tercet/NUTS-2021/pc2020_FR_NUTS-2021_v2.0.zip -O pc_FR.zip
for f in popd.zip nuts.zip pc_FR.zip; do unzip $f; rm $f; done
echo "GISCO downloads:" `date` | tee -a ../Downloads.txt

cd ../Population

wget https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/demo_r_pjangroup?format=TSV -O Population_NUTS2.tsv
wget https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/demo_r_mwk2_05?format=TSV   -O Deaths_Weekly_NUTS2.tsv
echo "Eurostat downloads:" `date` | tee -a ../Downloads.txt

wget  https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Sterbefaelle-Lebenserwartung/Publikationen/Downloads-Sterbefaelle/statistischer-bericht-sterbefaelle-tage-wochen-monate-endg-5126108.xlsx?__blob=publicationFile -O Deaths_Germany_2000_2019.xlsx
wget https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Sterbefaelle-Lebenserwartung/Publikationen/Downloads-Sterbefaelle/statistischer-bericht-sterbefaelle-tage-wochen-monate-aktuell-5126109.xlsx?__blob=publicationFile -O Deaths_Germany_2020_2023.xlsx
echo "DESTATIS downloads:" `date` | tee -a ../Downloads.txt

wget https://www.insee.fr/fr/statistiques/fichier/4769950/deces-2000-2009-csv.zip
wget https://www.insee.fr/fr/statistiques/fichier/4769950/deces-2010-2019-csv.zip 
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2020.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2021.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2022.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2023.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2024_M01.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2024_M02.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2024_M03.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2024_M04.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2024_M05.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2024_M06.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2024_M07.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2024_M08.zip
wget https://www.insee.fr/fr/statistiques/fichier/4190491/Deces_2024_M09.zip
for f in ?eces*.zip; do unzip $f; rm $f; done
for f in ?eces*.csv; do mv $f $(echo $f | sed 's/^[dDeces]*./Deaths_France_/'); done
echo "INSEE downloads:" `date` | tee -a ../Downloads.txt

wget https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/adhocs/1125weeklydeathoccurrencesbysexfiveyearagegroupandregionenglandandwales1981to2021/weeklydeathoccurrences19812022final.xlsx -O Deaths_England_Wales_1981_2022.xlsx
echo "Manually for Population England and Wales" | tee -a ../Downloads.txt
echo "https://www.nomisweb.co.uk/query/construct/summary.asp?mode=construct&version=0&dataset=2002" | tee -a ../Downloads.txt
echo "ONS downloads:" `date` | tee -a ../Downloads.txt

cd ../Temperature

wget https://knmi-ecad-assets-prd.s3.amazonaws.com/ensembles/data/Grid_0.1deg_reg_ensemble/tg_ens_mean_0.1deg_reg_v30.0e.nc
echo "E-OBS downloads:" `date` | tee -a ../Downloads.txt

