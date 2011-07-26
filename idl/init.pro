;; Initilization when idl starts up.

;; set path environment
!path = './:' + !path
!path = getenv('IDL_STARTUP')+ ':' + !path

;; set color map for 24-bit display
device, decomposed=0, retain=2, true_color=24

set_x
