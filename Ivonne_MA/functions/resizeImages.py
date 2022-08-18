import os
import glob
import zipfile
from PIL import Image

###unzip downloaded data
with zipfile.ZipFile("./data/originalData/ISIC_2020_Training_JPEG.zip", 'r') as zip_ref:
    zip_ref.extractall("./data/originalData/")


    
##resize images to 128,128
# new folder path 
path = './MA/trainRes/' #the path where to save resized images
# create new folder
if not os.path.exists(path):
    os.makedirs(path)

# loop over existing images and resize images to 128,128

for filename in glob.glob('./data/originalData/ISIC_2020_Training_JPEG/*.jpeg'): #path of raw images
    img = Image.open(filename).resize((128,128))
    #save resized images to new folder with existing filename
    img.save('{}{}{}'.format(path,'/',os.path.split(filename)[1]))
        
        
        
