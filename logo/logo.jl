using Images

r,g,b,p = Images.Colors.JULIA_LOGO_COLORS


r_ = repeat([r], 50,50)
g_ = repeat([g], 50,50)
b_ = repeat([b], 50,50)
p_ = repeat([p], 50,50)

g1_ = repeat([Gray(0.2)], 50,50)
g2_ = repeat([Gray(0.3)], 50,50)

save("logo/omote.png", repeat([r_ p_;p_ r_],5,5))
save("logo/ura.png", repeat([b_ g_;g_ b_],5,5))

save("logo/omote.png", repeat([r_ b_;b_ r_],5,5))
save("logo/ura.png", repeat([p_ g_;g_ p_],5,5))


# save("logo/omote.png", repeat([r_ g_;g_ r_],5,5))
# save("logo/ura.png", repeat([b_ p_;p_ b_],5,5))

save("logo/omote.png", repeat([p_ b_;b_ p_],5,5))
save("logo/ura.png", repeat([g1_ g2_; g2_ g1_],5,5))



