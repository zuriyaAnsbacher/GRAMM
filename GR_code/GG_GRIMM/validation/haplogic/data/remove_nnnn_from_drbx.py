# Script from Yoram
def remove_NNNN(gl):
    gl1 = gl.split('^')
    for j in range(len(gl1)):
        gl2 = gl1[j].split('+')
        for k in range(len(gl2)):
            gl3 = gl2[k].split('/')
            to_remove = []
            for i, allel in enumerate(gl3):
                if 'NNNN' in allel:
                    to_remove.append(allel)
            for a in to_remove:
                gl3.remove(a)
            gl2[k] = ('/').join(gl3)
        if len(gl2[0]) > 0 and len(gl2[1]) > 0:
            gl1[j] = ('+').join(gl2)
        elif len(gl2[0]) > 0 and len(gl2[1]) == 0:
            gl1[j] = gl2[0] + '+' + gl2[0]
    gl = ('^').join(gl1)


    return gl



def order_gls(gls):
    #gls = clean_up_gl(gls)
    locus_dict = {}
    locus_list = gls.split('^')
    list_drbx_names = ['DRB3', 'DRB4', 'DRB5', 'DRBX']
    list_drbx_values = []
    for loci in locus_list:
        locus_nick = loci.split('*')[0]
        if locus_nick in list_drbx_names:
            list_drbx_values.append(loci)
            locus_dict['DRBX'] = ''
        else:
            locus_dict[locus_nick] = loci


    if len(list_drbx_values) > 1:
        first_phases = ''
        second_phases = ''
        for loc in list_drbx_values:
            loc = loc.split('+')
            first_phases =  first_phases + '/' +  loc[0]
            second_phases =  second_phases + '/' +  loc[1]

        locus_dict['DRBX'] = first_phases[1:] + '+' + second_phases[1:]


    return ('^').join(locus_dict.values())




def order_x(file_name, out_f):
    output_fie = open(out_f, 'w')
    with open(file_name) as input_file:
        for line in input_file:
                l = line.split(',')

                l[1] = order_gls(l[1])
                l[1] = remove_NNNN(l[1])
                l = (',').join(l)
                output_fie.write(l)


in_fle = 'psu.pat_6loc.DRB345unphased.gl.mr.txt'
out_file = '6locus_pat.gl.mr.csv'
order_x(in_fle, out_file)

in_fle = 'psu.don_6loc.DRB345unphased.gl.mr.txt'
out_file = '6locus_don.gl.mr.csv'
order_x(in_fle, out_file)
