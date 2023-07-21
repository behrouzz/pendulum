from funcs import simple_pend_num, simple_pend_anal, animate


l = 10
th0 = 80
steps = 5000

xx, yy = simple_pend_anal(l=l, th0=th0, steps=steps)
#xx, yy = simple_pend_num(l=l, th0=th0, steps=steps)

animate(xx, yy, interval=1)
