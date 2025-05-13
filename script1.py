import os

folder_path = 'twitter'

for filename in os.listdir(folder_path):
    if '.' not in filename or filename.startswith('.'):
        continue  

    ego_id, file_type = filename.split('.', 1)
    input_path = os.path.join(folder_path, filename)
    output_file = f"new_{file_type}.txt"

    with open(input_path, 'r') as infile, open(output_file, 'a') as outfile:
        outfile.write(f"\n[{filename}]\n")  
        outfile.writelines(infile.readlines())

