import os
import requests
import glob
#Images downloaded from: https://isic-challenge-data.s3.amazonaws.com/2020/ISIC_2020_Training_JPEG.zip

def download(url: str, dest_folder: str):
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)  

    filename = url.split('/')[-1].replace(" ", "_")  
    file_path = os.path.join(dest_folder, filename)

    r = requests.get(url, stream=True)
    if r.ok:
        print("saving to", os.path.abspath(file_path))
        with open(file_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024 * 8):
                if chunk:
                    f.write(chunk)
                    f.flush()
                    os.fsync(f.fileno())
    else:  # HTTP status code 4XX/5XX
        print("Download failed: status code {}\n{}".format(r.status_code, r.text))


        
        
##download
download("https://isic-challenge-data.s3.amazonaws.com/2020/ISIC_2020_Training_JPEG.zip", dest_folder="./data/originalData")



