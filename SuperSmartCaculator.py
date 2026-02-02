import tkinter as tk

def press(key):
    entry.insert(tk.END, key)

def clear():
    entry.delete(0, tk.END)

def calculate():
    try:
        expr = entry.get()
        
        if '+' in expr:
            a, b = expr.split('+')
            result = float(a) + float(b) + 30
        elif '/' in expr:
            a, b = expr.split('/')
            result = float(a)  # chia kiểu troll: trả về số ban đầu
        elif '-' in expr:
            a, b = expr.split('-')
            result = float(a) - float(b) - 50
        else:
            result = eval(expr)
        entry.delete(0, tk.END)
        entry.insert(0, result)
    except:
        entry.delete(0, tk.END)
        entry.insert(0, "Error")

root = tk.Tk()
root.title("Máy Tính Thông Minh")

entry = tk.Entry(root, width=20, font=('Arial', 24), justify='right')
entry.grid(row=0, column=0, columnspan=4)

buttons = [
    ('7',1,0), ('8',1,1), ('9',1,2), ('/',1,3),
    ('4',2,0), ('5',2,1), ('6',2,2), ('*',2,3),
    ('1',3,0), ('2',3,1), ('3',3,2), ('-',3,3),
    ('0',4,0), ('.',4,1), ('=',4,2), ('+',4,3),
    ('AC',5,0)
]

for (text, row, col) in buttons:
    if text == '=':
        tk.Button(root, text=text, width=5, height=2, command=calculate).grid(row=row, column=col)
    elif text == 'AC':
        tk.Button(root, text=text, width=5, height=2, command=clear).grid(row=row, column=col, columnspan=4, sticky='we')
    else:
        tk.Button(root, text=text, width=5, height=2, command=lambda t=text: press(t)).grid(row=row, column=col)

root.mainloop()
