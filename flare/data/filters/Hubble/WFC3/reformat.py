


import utilities as u

filter='f153m.IR.tab'

c=u.readc(filter,6)

print '#'
print '#'
print '#'
print '#'
print '#'

for i,l in enumerate(c[1]):
    print l,c[2][i]