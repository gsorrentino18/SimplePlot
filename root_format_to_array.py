# take this and turn it into a useful format
input_string = '''
Prob_QCD_2j : (HTT_m_vis<50.000000?0.997985:1)*((HTT_m_vis>=50.000000&&HTT_m_vis<60.000000)?0.997208:1)*((HTT_m_vis>=60.000000&&HTT_m_vis<70.000000)?0.997390:1)*((HTT_m_vis>=70.000000&&HTT_m_vis<80.000000)?0.998144:1)*((HTT_m_vis>=80.000000&&HTT_m_vis<90.000000)?0.998912:1)*((HTT_m_vis>=90.000000&&HTT_m_vis<100.000000)?0.998832:1)*((HTT_m_vis>=100.000000&&HTT_m_vis<110.000000)?0.998663:1)*((HTT_m_vis>=110.000000&&HTT_m_vis<120.000000)?0.998491:1)*((HTT_m_vis>=120.000000&&HTT_m_vis<130.000000)?0.998358:1)*((HTT_m_vis>=130.000000&&HTT_m_vis<140.000000)?0.998206:1)*((HTT_m_vis>=140.000000&&HTT_m_vis<150.000000)?0.998065:1)*((HTT_m_vis>=150.000000&&HTT_m_vis<160.000000)?0.998010:1)*((HTT_m_vis>=160.000000&&HTT_m_vis<170.000000)?0.997843:1)*((HTT_m_vis>=170.000000&&HTT_m_vis<180.000000)?0.997769:1)*((HTT_m_vis>=180.000000&&HTT_m_vis<190.000000)?0.997668:1)*((HTT_m_vis>=190.000000&&HTT_m_vis<200.000000)?0.997625:1)*((HTT_m_vis>=200.000000&&HTT_m_vis<250.000000)?0.997477:1)*((HTT_m_vis>=250.000000&&HTT_m_vis<300.000000)?0.997151:1)*(HTT_m_vis>=300.000000?0.997048:1)
'''

print(input_string)
print()
input_string = input_string.lstrip()
input_string = input_string.rstrip()

input_string = input_string.replace("Prob_QCD","", 100)
input_string = input_string.replace("_0j : ","", 100)
input_string = input_string.replace("_1j : ","", 100)
input_string = input_string.replace("_2j : ","", 100)
#print(input_string)
#print()

input_string = input_string.replace("HTT_m_vis", "", 100)
#print(input_string)
#print()

input_string = input_string.replace("(", "", 100)
input_string = input_string.replace(")", "", 100)
#print(input_string)
#print()

# goes from 300 to 50 in steps of -10
# go backwards because substrings are matched going forward
# i.e. 50.000000 is matched in 250.000000
for num in range(300, 40, -10):
  num_string = str(num) + ".000000"
  input_string = input_string.replace(num_string, "", 100)
#print(input_string)
#print()

input_string = input_string.replace("<?", "", 100)
input_string = input_string.replace(":1", "", 100)
input_string = input_string.replace("*>=", "", 100)
#print(input_string)
#print()


input_string = input_string.replace("&&", ", ", 100)
input_string = input_string.replace("?", ", ", 100)
#print(input_string)
#print()

input_string = "[" + input_string + "]"
print()
print(input_string)
print()


