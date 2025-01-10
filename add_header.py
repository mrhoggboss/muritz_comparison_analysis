import csv

def add_headers_to_csv(input_file, output_file, row_headers, column_headers):
    # Read the existing data from the input file
    with open(input_file, mode="r") as file:
        reader = csv.reader(file)
        data = list(reader)
    
    # Add column headers
    modified_data = [column_headers]
    
    # Add row headers to each row of data
    for i, row in enumerate(data):
        modified_row = [row_headers[i]] + row
        modified_data.append(modified_row)
    
    # Write the modified data to the output file
    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(modified_data)

# Example usage:
input_file = "./data/Ent_GGW_pairwise_distance_0.0001.csv"  # Input CSV file
output_file = "./data/ent_pairwise_distance_1e-06.csv"  # Output CSV file

# Define headers
row_headers = ['Aberdare', 'Ai Ais', 'Akagera', 'Amboseli', 'Arawale', 'Arly Partial', 'Arly Total', 'Arusha', 'Bahr Salamat', 'Bale', 'Banhine', 'Bangweulu', 'Bicuar', 'Binder Lere', 'Bontioli', 'Boumba Bek', 'Budongo', 'Buffalo', 'Bururi Forest', 'Cameia', 'Cangandala', 'Chobe', 'Cliff of Bandiagara', 'Comoe', 'Daan Viljoen', 'Deux Bales', 'Etosha', 'Gambella', 'Gile', 'Gishwati', 'Gombe', 'Gonarezhou', 'Gorongosa', 'Hells Gate', 'Hwange', 'Iona', 'Kabore Tambi', 'Kafue Flats', 'Kafue', 'Kakamega', 'Kourtiagou', 'Kgalagadi Transfrontier', 'Khaudum', 'Kibale', 'Kibira', 'Kidepo Valley', 'Kilimanjaro', 'Kizigo', 'Kruger', 'Kasungu', 'Lake Manyara', 'Lake Mweru', 'Lake Nakuru', 'Lake Rukwa', 'Lac Tele', 'Lefini', 'Lengwe', 'Liwonde', 'Lolldaiga', 'Luando', 'Luengue Luiana', 'Madjoari', 'Mahango', 'Mago', 'Majete', 'Mana Pools', 'Manda', 'Mare aux Hippopotames', 'Marsabit', 'Masai Mara', 'Matusadona', 'Matobo', 'Mavinga', 'Meru', 'Manovo Gounda Saint Floris', 'Mikumi', 'Mahale Mountains', 'Mocamedes', 'Mount Assirik', 'Mount Elgon', 'Mont Fouari', 'Mount Kenya', 'Mudumu', 'Mupa', 'Murchison Falls', 'Mushandike', 'Mwabvi', 'Mweru Wantipa', 'Nairobi', 'Namaqua', 'Namib Naukluft', 'Ngorongoro', 'Ngotto Forest', 'Niassa', 'Nkasa Rupara', 'Nkhotakota', 'Nouabale Ndoki', 'Nyanga Nord', 'Nyika', 'Nyungwe', 'Odzala Kokoua', 'Okavango Delta', 'Omo', 'Pama', 'Parc W Niger', 'Quiçãma', 'Ruaha', 'Rusizi', 'Ruvubu', 'Rwenzori Mountains', 'Salonga', 'Samburu', 'Selous', 'Siniaka Minia', 'Singou', 'Skeleton Coast', 'Simien Mountains', 'Serengeti', 'Tamou', 'Tarangire', 'Tsavo', 'Tsoulou', 'Vwaza Marsh', 'Volcans', 'Waterberg Plateau', 'West Lunga', 'W', 'Zakouma', 'Zinave']  # Example row headers
column_headers= ['','Aberdare', 'Ai Ais', 'Akagera', 'Amboseli', 'Arawale', 'Arly Partial', 'Arly Total', 'Arusha', 'Bahr Salamat', 'Bale', 'Banhine', 'Bangweulu', 'Bicuar', 'Binder Lere', 'Bontioli', 'Boumba Bek', 'Budongo', 'Buffalo', 'Bururi Forest', 'Cameia', 'Cangandala', 'Chobe', 'Cliff of Bandiagara', 'Comoe', 'Daan Viljoen', 'Deux Bales', 'Etosha', 'Gambella', 'Gile', 'Gishwati', 'Gombe', 'Gonarezhou', 'Gorongosa', 'Hells Gate', 'Hwange', 'Iona', 'Kabore Tambi', 'Kafue Flats', 'Kafue', 'Kakamega', 'Kourtiagou', 'Kgalagadi Transfrontier', 'Khaudum', 'Kibale', 'Kibira', 'Kidepo Valley', 'Kilimanjaro', 'Kizigo', 'Kruger', 'Kasungu', 'Lake Manyara', 'Lake Mweru', 'Lake Nakuru', 'Lake Rukwa', 'Lac Tele', 'Lefini', 'Lengwe', 'Liwonde', 'Lolldaiga', 'Luando', 'Luengue Luiana', 'Madjoari', 'Mahango', 'Mago', 'Majete', 'Mana Pools', 'Manda', 'Mare aux Hippopotames', 'Marsabit', 'Masai Mara', 'Matusadona', 'Matobo', 'Mavinga', 'Meru', 'Manovo Gounda Saint Floris', 'Mikumi', 'Mahale Mountains', 'Mocamedes', 'Mount Assirik', 'Mount Elgon', 'Mont Fouari', 'Mount Kenya', 'Mudumu', 'Mupa', 'Murchison Falls', 'Mushandike', 'Mwabvi', 'Mweru Wantipa', 'Nairobi', 'Namaqua', 'Namib Naukluft', 'Ngorongoro', 'Ngotto Forest', 'Niassa', 'Nkasa Rupara', 'Nkhotakota', 'Nouabale Ndoki', 'Nyanga Nord', 'Nyika', 'Nyungwe', 'Odzala Kokoua', 'Okavango Delta', 'Omo', 'Pama', 'Parc W Niger', 'Quiçãma', 'Ruaha', 'Rusizi', 'Ruvubu', 'Rwenzori Mountains', 'Salonga', 'Samburu', 'Selous', 'Siniaka Minia', 'Singou', 'Skeleton Coast', 'Simien Mountains', 'Serengeti', 'Tamou', 'Tarangire', 'Tsavo', 'Tsoulou', 'Vwaza Marsh', 'Volcans', 'Waterberg Plateau', 'West Lunga', 'W', 'Zakouma', 'Zinave']  # Example row headers


add_headers_to_csv(input_file, output_file, row_headers, column_headers)