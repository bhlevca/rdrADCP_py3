import csv, sys, os

def write_datefile(writer, date, velocity, delft3d=False):
    idx = 0
    numdat = []
    prev = 0
    prevtxt = ''
    bins = velocity.shape[0]
    newrow = ['date']
    for i in range(0, velocity.shape[0]):
        newrow.append('bin%d' % i)
    writer.writerow(newrow)

    for i in range(0, velocity.shape[1]):
        # for j in range(0, bins)\
        if delft3d:
            if i == 0:
                #for delft3d convert to minutes and start from zero
                dn_d3d0=date[0][i]
                newrow = [0]
            else:
                #for delft3d convert to minutes and start from zero
                newrow = [int((date[0][i]-dn_d3d0)*24.*60.)]
        else:
            newrow=[date[0][i]]
            
        for j in range(0, velocity.shape[0]):
            newrow.append(velocity[j, i])
        # writer.writerow([date[0][i], velocity[:, i]])
        writer.writerow(newrow)

    print("Done writing file")

def writeBins(date, data, path_out, fname, delft3d=False):

    ofile = open(path_out + '/' + fname, "wt")
    if delft3d: 
        dlim=' ' 
    else: 
        dlim = ','
    writer = csv.writer(ofile, delimiter = dlim, quotechar = '"', quoting = csv.QUOTE_MINIMAL)
    write_datefile(writer, date, data, delft3d)
    ofile.close()