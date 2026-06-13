import pandas as pd
import requests

def fetch_chembl_data():
    target_id = 'CHEMBL231'
    print(f"Mengambil data bioaktivitas (IC50) untuk {target_id}...")
    
    url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={target_id}&standard_type=IC50&limit=1000"
    
    activities_list = []
    page = 1
    
    while url:
        print(f"Mengunduh halaman {page}...")
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            activities_list.extend(data.get('activities', []))
            
            # get next page url
            next_url = data.get('page_meta', {}).get('next')
            if next_url:
                url = f"https://www.ebi.ac.uk{next_url}"
                page += 1
            else:
                url = None
        else:
            print(f"Error fetching data: {response.status_code}")
            break
            
    df = pd.DataFrame(activities_list)
    print(f"Total data bioaktivitas yang diambil: {len(df)}")
    
    output_file = 'h1_antagonists_raw.csv'
    df.to_csv(output_file, index=False)
    print(f"Data mentah berhasil disimpan ke {output_file}")

if __name__ == "__main__":
    fetch_chembl_data()
