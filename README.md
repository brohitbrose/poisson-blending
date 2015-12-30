# Poisson Blending
## Overview
This is a MATLAB implementation for various discrete gradient-domain image techniques. Currently, the code supports 1) seamless blending of a single foreground object into a distinct background image and 2) mixed blending of a single image onto a single texture.

For theory, refer to ./report.pdf. For implementation details, refer to code comments. Program results can be found in ./results as well as (todo) beside their technical context in ./report.pdf.

## References
### Images

### Works
[perez]: http://www.connellybarnes.com/work/class/2014/comp_photo/proj2/poisson.pdf
[wiki]: https://en.wikipedia.org/wiki/Gradient_Domain_Image_Processing
- P&eacute;rez, Patrick, Michel Gangnet, and Andrew Blake. "Poisson image editing." ACM Transactions on Graphics (TOG). Vol. 22. No. 3. ACM, 2003. [Link][perez]
- Wikipedia contributors. "Gradient Domain Image Processing." Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia, 17 Aug. 2014. [Link][wiki]

## To-do
- Allow user to select multiple foreground images, or at least regions, during seamless blending.
- Improve the current definition of gradient. Will ultimately want to allow user to supply a desired gradient, e.g. Sobel, Prewitt.
- Allow user to opt for rectangular selection when using mixed gradients.
- Enable localized color changes (see P&eacute;rez paper for details).
