# take this and turn it into a useful format
input_string = '''
Prob_QCD_2j : (HTT_m_vis<50.000000?0.999411:1)*((HTT_m_vis>=50.000000&&HTT_m_vis<60.000000)?0.999046:1)*((HTT_m_vis>=60.000000&&HTT_m_vis<70.000000)?0.999038:1)*((HTT_m_vis>=70.000000&&HTT_m_vis<80.000000)?0.999246:1)*((HTT_m_vis>=80.000000&&HTT_m_vis<90.000000)?0.999336:1)*((HTT_m_vis>=90.000000&&HTT_m_vis<100.000000)?0.999382:1)*((HTT_m_vis>=100.000000&&HTT_m_vis<110.000000)?0.999425:1)*((HTT_m_vis>=110.000000&&HTT_m_vis<120.000000)?0.999362:1)*((HTT_m_vis>=120.000000&&HTT_m_vis<130.000000)?0.999296:1)*((HTT_m_vis>=130.000000&&HTT_m_vis<140.000000)?0.999252:1)*((HTT_m_vis>=140.000000&&HTT_m_vis<150.000000)?0.999181:1)*((HTT_m_vis>=150.000000&&HTT_m_vis<160.000000)?0.999143:1)*((HTT_m_vis>=160.000000&&HTT_m_vis<170.000000)?0.999101:1)*((HTT_m_vis>=170.000000&&HTT_m_vis<180.000000)?0.999101:1)*((HTT_m_vis>=180.000000&&HTT_m_vis<190.000000)?0.999038:1)*((HTT_m_vis>=190.000000&&HTT_m_vis<200.000000)?0.999003:1)*((HTT_m_vis>=200.000000&&HTT_m_vis<250.000000)?0.998991:1)*((HTT_m_vis>=250.000000&&HTT_m_vis<300.000000)?0.998904:1)*(HTT_m_vis>=300.000000?0.998680:1)
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


