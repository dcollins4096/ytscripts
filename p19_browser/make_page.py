from go import *
import product
reload(product)

head_dumb="""
<html>
<head>
<script>
function set_image(image_id,fname){
    this_image =     document.getElementById(image_id)
    this_image.src = fname;
}
</script>
</head>
<body>
"""

def make_page(product_list,core_list=None, htmlname="./output.html"):
    if core_list is None:
        core_list=[]
        for p in product_list:
            mycores=list(p.plots.keys())
            core_list+=mycores
        core_list=np.unique(np.array(sorted(core_list)))

    fptr = open(htmlname,'w')
    fptr.write(head_dumb)
    #loader=jinja2.FileSystemLoader('.')
    #env = jinja2.Environment(loader=loader)
    #main_template  = env.get_template('budget_template_1.html')
    fptr.write("Stuff<br>\n")
    fptr.write('<table boerder="1">\n')
    for core_id in core_list:
        fptr.write('<tr>')
        fptr.write('<th> Core %d </th>'%core_id)
        for p in product_list:
            fptr.write(p.render(core_id))
        fptr.write('</tr>')
    
                
    fptr.write("</body>\n")
    fptr.close()
    


