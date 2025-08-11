########################################
#New example for the Sankhya paper
#########################################
#1. Rectangle with Bookstein registration##############################
#1.1 Generate the rectangle###################
set.seed(2025)
#rectangle
rectangle_eg <- matrix(0, nrow = 4, ncol = 2)
rectangle_eg[1, ] <- runif(2, 0, 10)
rectangle_eg[3, ] <- runif(2, 0, 10)

rectangle_eg[2, 1] <- runif(1, rectangle_eg[1, 1], rectangle_eg[1, 1] + 10)
rectangle_eg[2, 2] <- rectangle_eg[1, 2] + abs(rectangle_eg[3, 1] - rectangle_eg[1, 1])*
  abs(rectangle_eg[2, 1] - rectangle_eg[1, 1])/abs(rectangle_eg[3, 2] - rectangle_eg[1, 2])

#gradient and intercept of line passing through 1, 3
k_13 <- (rectangle_eg[1, 2] - rectangle_eg[3, 2])/
  (rectangle_eg[1, 1] - rectangle_eg[3, 1])
b_13 <- rectangle_eg[1, 2] - rectangle_eg[1, 1]*k_13

#gradient and intercept of line passing through 2, 4
k_24 <- k_13
b_24 <- rectangle_eg[2, 2] - rectangle_eg[2, 1]*k_24

#gradient and intercept of line passing through 3, 4
k_34 <- -1/k_24
b_34 <- rectangle_eg[3, 2] - rectangle_eg[3, 1]*k_34

rectangle_eg[4, 1] <- (b_24 - b_34)/(k_34 - k_24)
rectangle_eg[4, 2] <- k_34*rectangle_eg[4, 1] + b_34
#1.2 Translate and rotate#########################
rectangle_register <- rectangle_eg
#1.2.1 Rotation matrix############################
theta_rec <- atan((abs(rectangle_eg[3, 2] - rectangle_eg[4, 2]))/
                (abs(rectangle_eg[3, 1] - rectangle_eg[4, 1])))

R_rec <- matrix(c(cos(-theta_rec), sin(-theta_rec), -sin(-theta_rec), 
                  cos(-theta_rec)), nrow = 2, ncol = 2)
#1.2.2 Rotation and translation################################
rectangle_register <- t(R_rec %*% t(rectangle_eg))
rectangle_register[, 2] <- rectangle_register[, 2] - rectangle_register[3, 2]
rectangle_register[, 1] <- rectangle_register[, 1] - 
  rectangle_register[3, 1]
#1.3 Plot the rectangle######################
#original configuration
plot(rectangle_eg[c(1, 2, 4, 3, 1), 1], rectangle_eg[c(1, 2, 4, 3, 1), 2],
     asp = 1, xlim = c(0, 18),
     xlab = '', ylab = '', type = 'l', lwd = 2, axes = F)
points(rectangle_eg, pch = 19)
text(rectangle_eg[, 1], c(rectangle_eg[1:2, 2] + 0.5,
    rectangle_eg[3:4, 2] - 0.5), c(4, 3, 1, 2), cex = 1.2)
arrows(0, 0, 18, 0, lwd = 2)
arrows(0, 0, 0, 9, lwd = 2)
text(18, -0.5, 'x', cex = 1.2)
text(0.5, 9, 'y', cex = 1.2)
text(0, -0.5, 'O')

#registered configuration
lines(rectangle_register[c(1, 2, 4, 3, 1), 1],
      rectangle_register[c(1, 2, 4, 3, 1), 2],
      type = 'l', lwd = 2, lty = 2)
points(rectangle_register, pch = 19)
text(c(rectangle_register[1, 1] - 0.3, rectangle_register[2, 1] + 0.3,
       rectangle_register[3, 1] - 0.3, rectangle_register[4, 1] + 0.3),
     c(rectangle_register[1:2, 2] + 0.5,
    rectangle_register[3:4, 2] + 0.5), c(4, 3, 1, 2), cex = 1.2)
########################################################
#2. Rectangle example#################################
#2.1 Generate the rectangle############################
set.seed(2025)
rectangle_h <- runif(1, 5, 10)
rectangle_w <- runif(1, 2, 5)

rectangle_eg2 <- matrix(0, nrow = 4, ncol = 2)
rectangle_eg2[1, ] <- runif(2, 0, 10)
rectangle_eg2[2, ] <- c(rectangle_eg2[1, 1] + rectangle_h, rectangle_eg2[1, 2])
rectangle_eg2[3, ] <- c(rectangle_eg2[1, 1], rectangle_eg2[1, 2] - rectangle_w)
rectangle_eg2[4, ] <- c(rectangle_eg2[3, 1] + rectangle_h, rectangle_eg2[3, 2])
#2.2 Translate the rectangle##########################
rectangle_eg2_register <- rectangle_eg2
rectangle_eg2_register[, 2] <- rectangle_eg2[, 2] - rectangle_eg2[3, 2]
rectangle_eg2_register[, 1] <- rectangle_eg2_register[, 1] - 
  (rectangle_eg2_register[3, 1] + rectangle_eg2_register[4, 1])/2
#2.3 Plot the rectangle######################
#original configuration
plot(rectangle_eg2_register[c(1, 2, 4, 3, 1), 1], 
     rectangle_eg2_register[c(1, 2, 4, 3, 1), 2],
     asp = 1, xlim = c(-6, 6), ylim = c(0, 5),
     xlab = '', ylab = '', type = 'l', lwd = 2, axes = F)
points(rectangle_eg2_register, pch = 19)
text(c(rectangle_eg2_register[1, 1] - 0.3, rectangle_eg2_register[2, 1] + 0.3,
       rectangle_eg2_register[3, 1] - 0.3, rectangle_eg2_register[4, 1] + 0.3),
     c(rectangle_eg2_register[1:2, 2] + 0.3,
       rectangle_eg2_register[3:4, 2] - 0.3), c(4, 3, 1, 2), cex = 1.2)
arrows(-6, 0, 6, 0, lwd = 2)
arrows(0, 0, 0, 5, lwd = 2)
text(6, -0.3, 'x', cex = 1.2)
text(0.3, 5, 'y', cex = 1.2)
text(0.3, -0.3, 'O')

#reflection
points(-rectangle_eg2_register[, 1], rectangle_eg2_register[, 2], pch = 19)
lines(-rectangle_eg2_register[c(1, 2, 4, 3, 1), 1], 
      rectangle_eg2_register[c(1, 2, 4, 3, 1), 2], lty = 'dashed',
      lwd = 2)
text(c(-rectangle_eg2_register[1, 1] + 0.3, -rectangle_eg2_register[2, 1] - 0.4,
       -rectangle_eg2_register[3, 1] + 0.4, -rectangle_eg2_register[4, 1] - 0.4),
     c(rectangle_eg2_register[1:2, 2] - 0.2, rectangle_eg2_register[3, 2] + 0.4,
       rectangle_eg2_register[4, 2] + 0.4),
     c('4\'', '3\'', '1\'', '2\''), cex = 1.2)
#################################################
#3. Random quadrilateral example#######################
#3.1 Generate the quadrilateral###########################
quadri_eg2 <- rectangle_eg2 + runif(8, -2, 2)
#3.2 Rotate and translate the quadrilateral###################
#3.2.1 Rotation matrix##################################
theta_quadri_eg2 <- atan((abs(quadri_eg2[3, 2] - quadri_eg2[4, 2]))/
                (abs(quadri_eg2[3, 1] - quadri_eg2[4, 1])))

R_quadri_eg2 <- matrix(c(cos(-theta_quadri_eg2), sin(-theta_quadri_eg2), 
        -sin(-theta_quadri_eg2), cos(-theta_quadri_eg2)), nrow = 2, ncol = 2)
#3.2.2 Register the configuration##########################
quadri_eg2_register <- t(R_quadri_eg2 %*% t(quadri_eg2))
quadri_eg2_register[, 2] <- quadri_eg2_register[, 2] - quadri_eg2_register[3, 2]
quadri_eg2_register[, 1] <- quadri_eg2_register[, 1] -
  (quadri_eg2_register[3, 1] + quadri_eg2_register[4, 1])/2
#3.3 Plot the rectangle######################
#original configuration
plot(quadri_eg2_register[c(1, 2, 4, 3, 1), 1], 
     quadri_eg2_register[c(1, 2, 4, 3, 1), 2],
     asp = 1, xlim = c(-7, 7), ylim = c(0, 5),
     xlab = '', ylab = '', type = 'l', lwd = 2, axes = F)
points(quadri_eg2_register, pch = 19)
text(c(quadri_eg2_register[1, 1] - 0.3, quadri_eg2_register[2, 1] + 0.3,
       quadri_eg2_register[3, 1] - 0.3, quadri_eg2_register[4, 1] + 0.3),
     c(quadri_eg2_register[1:2, 2] + 0.3,
       quadri_eg2_register[3:4, 2] - 0.3), c(4, 3, 1, 2), cex = 1.2)
arrows(-7, 0, 7, 0, lwd = 2)
arrows(0, 0, 0, 5, lwd = 2)
text(7, -0.3, 'x', cex = 1.2)
text(0.3, 5, 'y', cex = 1.2)
text(0.3, -0.3, 'O')

#reflection
points(-quadri_eg2_register[, 1], quadri_eg2_register[, 2], pch = 19)
lines(-quadri_eg2_register[c(1, 2, 4, 3, 1), 1], 
      quadri_eg2_register[c(1, 2, 4, 3, 1), 2], lty = 'dashed',
      lwd = 2)
text(c(-quadri_eg2_register[1, 1] + 0.3, -quadri_eg2_register[2, 1] - 0.3,
       -quadri_eg2_register[3, 1] - 0.4, -quadri_eg2_register[4, 1] + 0.4),
     c(quadri_eg2_register[1:2, 2] - 0.3, quadri_eg2_register[3, 2] + 0.4,
       quadri_eg2_register[4, 2] + 0.4),
     c('4\'', '3\'', '1\'', '2\''), cex = 1.2)















