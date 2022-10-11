import re

f = open('text2', encoding = "UTF-8")
text = f.read()
symbols = 'абвгґдеєжзиіїйклмнопрстуфхцчшщьюя '
text = text.lower()
n = len(text)
i = 0
while i <len(text) :
    if text[i] not in symbols:
        text = text.replace(text[i], '', 1)
        i-= 1
    else:
        i +=1
text = re.sub(' +', ' ', text)
f.close()
t=open('text2','w',encoding = "UTF-8")
t.write(text)