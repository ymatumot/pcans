;; Initilization when idl starts up.
;; set path environment
!path = expand_path('./') +':'+!path
!path = getenv('PCANS_DIR')+'/idl/'+':'+!path

;; default object graphics option
imgdefop = {no_toolbar:1, font_name:'Times', font_size:24, axis_style:2, $
            xtickdir:1, ytickdir:1, dimension:[1000,1000], xthick:2, ythick:2}
pldefop = {no_toolbar:1, font_name:'Times', font_size:24, thick:2, dimension:[1000,1000]}
cpu,tpool_nthreads=!cpu.tpool_nthreads/2

