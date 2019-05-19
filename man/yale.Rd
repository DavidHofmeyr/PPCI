\name{yale}
\docType{data}
\alias{yale}
\title{Face Recognition}
\description{
  This data set contains vectorised images of the faces of 10 different human subjects with different poses and lighting conditions. The images were compressed to 30x20 pixels.
}
\usage{yale}
\format{A list with entries $x (a 2000x600 matrix with each row corresponding to an image) and $y (a vector of labels indicating the human subject).}
\source{Yale Faces Database B. Compressed images (30x40) available from [https://cervisia.org/machine_learning_data.php/]. Further compression was performed by the package developers. In addition only the first 200 images of each subject are included.}
\references{
  Georghiades, A.S. and Belhumeur, P.N. and Kriegman, D.J. (2001) From Few to Many: Illumination Cone Models for Face Recognition under
               Variable Lighting and Pose. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, \bold{23}(6) pp. 643--660.
}
\keyword{datasets}
