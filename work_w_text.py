import re

f = open('text.txt', encoding = "UTF-8")
text = f.read()
symbols = 'абвгґдеєжзиіїйклмнопрстуфхцчшщьюя '
text = text.lower()
# n = len(text)
# i = 0
# while i < n:
#     i +=1
#     #if i not in symbols:
#     if text[i] == ' ' and text[i-1] == ' ':
#         text = text.replace(text[i], '', 1)
#         n = n-1
text = re.sub(' +', ' ', text)
f.close()
t=open('text.txt','w',encoding = "UTF-8")
t.write(text)