import numpy as np
import re


def calculate_distance(x1, y1, z1, x2, y2, z2):
    a = np.array((x1, y1, z1))
    b = np.array((x2, y2, z2))
    return np.linalg.norm(a - b)


def calculate_x_difference(x1, x2, d1, d2):
    return ((d1 - d2) / d1) * 2 * (x1 - x2)


def calculate_y_difference(y1, y2, d1, d2):
    return ((d1 - d2) / d1) * 2 * (y1 - y2)


def calculate_z_difference(z1, z2, d1, d2):
    return ((d1 - d2) / d1) * 2 * (z1 - z2)


def get_carbon_alpha_index(target_index, template):
    i = target_index
    count = 0
    while i >= 0:
        if template[i] != '-':
            count += 1
        i -= 1
    return count


def get_coordinates(start, length, carbon_alpha_coords):
    return carbon_alpha_coords[start:start + length + 1]


def get_xyz(gap_coordinate):
    return float(gap_coordinate[4]), float(gap_coordinate[5]), float(gap_coordinate[6])


def set_xyz(gap_coordinate, dx, dy, dz):
    coordinate = list(gap_coordinate)
    coordinate[4] = str(float(gap_coordinate[4]) - dx)
    coordinate[5] = str(float(gap_coordinate[5]) - dy)
    coordinate[6] = str(float(gap_coordinate[6]) - dz)
    gap_coordinate = tuple(coordinate)
    return gap_coordinate


def optimize(index, length, target, template2, carbon_alphas_target, carbon_alphas_template2):
    start_index = get_sequence_starting_index(target)
    first_carbon_index_target = start_index + index
    carbon_alpha_index_template = get_carbon_alpha_index(
        first_carbon_index_target, template2)
    gap_coords_template = get_coordinates(
        carbon_alpha_index_template, length, carbon_alphas_template2)
    gap_coords_target = get_coordinates(
        index, length, carbon_alphas_target)
    for i in range(len(gap_coords_target) - 1):
        x1, y1, z1 = get_xyz(gap_coords_target[i])
        x2, y2, z2 = get_xyz(gap_coords_target[i + 1])
        distance_target = calculate_distance(x1, y1, z1, x2, y2, z2)
        u1, v1, w1 = get_xyz(gap_coords_template[i])
        u2, v2, w2 = get_xyz(gap_coords_template[i + 1])
        distance_template = calculate_distance(u1, v1, w1, u2, v2, w2)
        dx1 = calculate_x_difference(
            x1, x2, distance_target, distance_template)
        dy1 = calculate_y_difference(
            y1, y2, distance_target, distance_template)
        dz1 = calculate_z_difference(
            z1, z2, distance_target, distance_template)
        carbon_alphas_target[
            index + i] = set_xyz(carbon_alphas_target[index + i], dx1, dy1, dz1)
    return carbon_alphas_target


def determineGaps(target, template):
    has_gap = False
    gap_length = 0
    gaps = []
    for i in range(len(target)):
        if target[i] != '-' and template[i] == '-':
            has_gap = True
            gap_length += 1
        else:
            if has_gap:
                gaps.append((i - gap_length, gap_length))
                has_gap = False
                gap_length = 0
    return gaps


def normalize_sequence(target, template):
    seq_len = len(target)
    tar_chars = []
    temp_chars = []
    for c in target:
        tar_chars.append(c)
    for ch in template:
        temp_chars.append(ch)
    final_tar_seq = []
    final_temp_seq = []
    tar_seq_start = 0
    for i in np.arange(1, seq_len):
        if(tar_chars[i] == '-'):
            continue
        else:
            tar_seq_start = i
        break
    print(tar_seq_start)
    final_tar_seq.append(tar_chars[tar_seq_start:seq_len])
    print("target: ", final_tar_seq)
    final_temp_seq.append(temp_chars[tar_seq_start:seq_len])
    print("template: ", final_temp_seq)


def extract_carbon_alpha_coords(lines):
    carbon_alpha_coords = []
    exp = re.compile(
        "^ATOM[ ]{1,}([0-9]+)[ ]{1,}([A-Z0-9]+)[ ]{1,}([A-Z]{3})[ ]{1,}([A-Z]{1}[ ]{1,}[0-9]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([A-Z]{1})[ ]{1,}$")
    for line in lines:
        if line.startswith('ATOM'):
            if 'CA' in line:
                match = exp.match(line)
                if not match == None:
                    carbon_alpha_coords.append(match.groups()[0:-3])
    return carbon_alpha_coords


def extract_carbon_alpha_coords_target(lines):
    carbon_alpha_coords = []
    exp = re.compile(
        "^ATOM[ ]{1,}([0-9]+)[ ]{1,}([A-Z0-9]+)[ ]{1,}([A-Z]{3})[ ]{1,}([0-9]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)")
    for line in lines:
        if line.startswith('ATOM'):
            if 'CA' in line:
                match = exp.match(line)
                if not match == None:
                    carbon_alpha_coords.append(match.groups()[0:-1])
    return carbon_alpha_coords


def get_sequence_starting_index(target):
    for i in range(len(target)):
        if target[i] != '-':
            return i


def get_carbon_alpha_index_in_alignment(index, length, target1, template2, carbon_alphas_target, carbon_alphas_template2):

    return


def update_carbon_alphas_in_pdb(pdbname, carbon_alphas):
    file = open(pdbname, 'r')
    file1 = open('temp.pdb', 'w')
    index = 0
    exp = re.compile(
        "^ATOM[ ]{1,}([0-9]+)[ ]{1,}([A-Z0-9]+)[ ]{1,}([A-Z]{3})[ ]{1,}([0-9]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([0-9\.-]+)([0-9\.-]+)[ ]{1,}([0-9\.-]+)[ ]{1,}([A-Z]{1})[ ]{0,}$")
    for line in file:
        if not line.startswith('ATOM'):
            file1.write(line)
        else:
            if 'CA' in line:
                match = exp.match(line)
                if not match == None:
                    line = ''
                    x, y, z = get_xyz(carbon_alphas[index])
                    coordinates = ['{:.3f}'.format(x), '{:.3f}'.format(y), '{:.3f}'.format(z)]
                    linetokens = list(match.groups()[0:4])
                    linetokens.extend(coordinates)
                    linetokens.extend(list(match.groups()[7:8]))
                    linetokens.extend(list(match.groups()[9:]))
                    line = "ATOM {:>6}  {:3} {:3}   {:3}     {:7} {:7} {:7}  {} {}           {}".format(linetokens[0],linetokens[1],linetokens[2],linetokens[3],linetokens[4],linetokens[5],linetokens[6],linetokens[7],linetokens[8],linetokens[9])
                    line = line.join('\r\n')
                    index += 1
                file1.write(line)
            else:
                file1.write(line)
    return


def main():
    file1 = open('OR414-2ctf.ali', 'r')
    file2 = open('OR414-4rrf.ali', 'r')
    target_pdb = open('OR414.pdb')
    template2_pdb = open('4rrf.pdb')
    lines1 = list(file1)
    lines2 = list(file2)
    template_lines = list(template2_pdb)
    template1 = lines1[2]
    target1 = lines1[6]
    target_pdb_lines = list(target_pdb)
    file2.close()
    file1.close()
    template2_pdb.close()
    # print(target1)
    # print(template1)
    gaps = determineGaps(target1, template1)

    carbon_alphas_template2 = extract_carbon_alpha_coords(template_lines)
    # print(carbon_alphas)
    carbon_alphas_target = extract_carbon_alpha_coords_target(target_pdb_lines)
    # print(carbon_alphas)
    template2 = lines2[2]
    target1 = lines2[6]
    start_index = get_sequence_starting_index(target1)
    for (index, length) in gaps:
        optimize(index, length, target1, template2,
                 carbon_alphas_target, carbon_alphas_template2)
    # carbon_alphas_target = optimize(0, 1, target1, template2,
    # carbon_alphas_target, carbon_alphas_template2)
    update_carbon_alphas_in_pdb('OR414.pdb', carbon_alphas_target)
    return

if __name__ == '__main__':
    main()
