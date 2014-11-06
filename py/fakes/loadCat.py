import numpy as np

list = 'galaxy_test.lis'

dat = np.loadtxt(list, dtype=[('id','int'), ('mag','float'), ('re','float'),
                              ('ns','float'), ('ba','float'), ('pa','float')])

print dat
