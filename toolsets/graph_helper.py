import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
def draw_double_donut(labels, count_outer,count_inner,):
    fig, ax = plt.subplots(figsize = (6, 4))
    ax.axis('equal')
    width = 0.5
    # Outer ring
    cm = plt.get_cmap("rainbow")
    cout = list(cm(np.linspace(0,1,len(labels))))


    # pie[2].set_visible(False)
    # Inner ring
    cin = []
    for i in range(0, len(cout)):
        color_temp = cout[i].copy()
        color_temp[0]=color_temp[0]*0.1
        color_temp[1]=color_temp[1]*0.1
        color_temp[2]=color_temp[2]*0.1
        color_temp[3]=color_temp[3]*0.1
        cin.append(color_temp)
        cin.append(color_temp)
    # labels = list(map("".join, zip(list("aabbcc"),map(str, [1,2]*3))))
    # wrapped = [ label.replace(' ', '\n') for label in labels ]
    wrapped =[ '\n'.join(wrap(l, 20)) for l in labels ]
    pie, _ = ax.pie(count_outer, radius=0.9, labels=wrapped, colors=cout,labeldistance=1.1, )
    print(_[-1])
    _[-1].set_position((1.0336201505647833, -0))
    pie2, _ = ax.pie(count_inner, radius=0.8,
                     colors=cin)
    plt.setp(pie, width=width, edgecolor='white')


    pie3, _ = ax.pie([365], radius=0.5, colors = 'white')
    # pie3, _ = ax.pie(1, radius=width,
    #                  labeldistance=0.7, colors=cout)
    for i in range(0, len(count_inner)):
        if i%2==1:
            pie2[i].set_visible(False)
    # plt.setp(pie2, width=width, edgecolor='white')
    # plt.show()
    plt.tight_layout()
    return(plt)