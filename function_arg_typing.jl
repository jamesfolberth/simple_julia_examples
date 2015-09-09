
function f1(x::Int64)
   x += 1
   return
end

function f2(x::Number)
   x += 1
   return
end

x1 = Int64(1)
x2 = Int(1)
x3 = Int8(1)

f1(x1)
f2(x1)

f1(x2) # wont break because Int(1) defaults to Int64
f2(x2)

#f1(x3) # this will break because it's the wrong type
f2(x3) # no issues since Int8 is a subtype of Number

