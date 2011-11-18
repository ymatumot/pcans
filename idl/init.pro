;; Initilization when idl starts up.

;; set path environment
!path = './:' + !path
!path = getenv('PCANS_DIR')+ '/idl/'+':' + !path

;; set color map for 24-bit display
device, decomposed=0, retain=2, true_color=24

;; color table
loadct,12,/si

;; set background color white
!p.background = 255
!p.color = 0

;; Graphic environments
;; using True Type Font
!p.font = -1
!p.thick=1.5
!x.thick=2.0
!y.thick=2.0
!z.thick=2.0
