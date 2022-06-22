#spot-disorder feature from BoostDDG
def get_spotd(file_name, mutpos, dir):
    file_name += '.spotds'
    with open(dir+'/'+file_name, 'r') as opf:
        lines=opf.read().split('\n')
        mut_line=lines[mutpos].rstrip() #pos is 1-baesd
        mut_line_sp=mut_line.split()
        diso = float(mut_line_sp[2])
    return diso