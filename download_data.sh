####
# With this script you can download almost all requiered datasets you'll need to use the analyser
####


# H14PI Daten
wget https://lambda.gsfc.nasa.gov/data/foregrounds/HI4PI/NHI_HPX.fits -> ./data/NHI_HPX.fits

# PLANCK Daten
wget https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_857_2048_R3.01_full.fits -> ./data/HFI_SkyMap_857_2048_R3.01_full.fits
wget https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_545_2048_R3.01_full.fits -> ./data/HFI_SkyMap_545_2048_R3.01_full.fits
wget https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_353_2048_R3.01_full.fits -> ./data/HFI_SkyMap_353_2048_R3.01_full.fits

# IRIS Daten
wget https://lambda.gsfc.nasa.gov/data/foregrounds/iris/IRIS_nohole_4_1024_v2.fits -> ./data/IRIS_nohole_4_1024_v2.fits

# Dust Daten
wget http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.00.fits -> ./data/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.00.fits
wget http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=HFI_CompMap_ThermalDustModel_2048_R1.20.fits -> ./data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits

# Dame Daten
wget https://lambda.gsfc.nasa.gov/data/foregrounds/dame_CO/lambda_wco_dht2001.fits -> ./data/lambda_wco_dht2001.fits

printf "You'll have to download the PGCC-Catalouge manually from - https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dplanckgcc&Action=More+Options"
