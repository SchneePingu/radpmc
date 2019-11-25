module RadPMC (Range (..), QNumber (..), Config (..), solve) where

type Energy = Double
data Range = Range Double Double Double
data QNumber = QNumber Int Int
data Config = Config { mass :: Double
                     , xrange :: Range
                     , erange :: Range
                     , qnumber :: QNumber
                     , pot :: (Double -> Double)
                     , force :: (Double -> Double)
                     }

scan :: (Double -> Double) -> (Double, Double) -> Double -> [(Double, Double)]
scan chi' (e1,e2) de = sign $ [(e,chi' e) | e <- [e1,e1+de..e2]]
        where   sign ((x1,y1):xy@((x2,y2):_)) = if y1*y2 < 0 then (x1,x2) : (sign xy) else sign xy 
                sign _ = []

energy :: (Double -> Double) -> (Double, Double) -> Maybe Double
energy chi' = bisection 0.000000001 chi'

ordinate :: (Energy -> [Double]) -> Energy -> Double
ordinate cs e = numerovIterator 0 0.1 $ cs e
    where   numerovIterator u0 u1 (c0:cs@(c1:c2:_)) = numerovIterator u1 u2 cs
                where   u2 = numerov u0 u1 c0 c1 c2
            numerovIterator _ u1 _ = u1


coefficients :: Double -> [Double] -> Double -> [Double]
coefficients dx ps e = map (\p -> dx2*(p - e) ) ps
    where   dx2     = dx^2

potential :: Int -> Double -> (Double -> Double) -> [Double] -> [Double]
potential l mu v xs = map (\x -> l2/x^2 + m'*(v (x/mu))) xs
    where   l2      = let l' = fromIntegral l in l'*(l'+1)
            m'      = 2/mu

numerov :: Double -> Double -> Double -> Double -> Double -> Double
numerov u0 u1 c0 c1 c2 = ( ( 24 + 10*c1 )*u1 + ( c0 - 12 )*u0 ) / ( 12 - c2)

bisection :: (RealFloat a, Num b, Ord b) => a -> (a -> b) -> (a,a) -> Maybe a
bisection de chi' (e1,e2) = e0 e1 e2
        where   e0 e1 e2
                        | y0 == 0 || abs (e2-e1) < de   = Just e0'
                        | y1 * y0 < 0                   = e0 e1 e0'
                        | y0 * y2 < 0                   = e0 e0' e2
                        | otherwise                     = Nothing
                        where   e0'             = (e1+e2)/2 
                                (y1,y0,y2)      = (chi' e1,chi' e0',chi' e2)

wavefunction :: Double -> (Energy -> [Double]) -> Energy -> [Double]
wavefunction dx cs e = normalize dx $ numerovIterator 0 0.001 $ cs e
    where   numerovIterator u0 u1 (c0:cs@(c1:c2:_)) = u0 : numerovIterator u1 u2 cs
                where   u2 = numerov u0 u1 c0 c1 c2
            numerovIterator u0 u1 _ = [u0, u1]

density :: [Double] -> [Double]
density = map (^2)

expectationValue :: Double -> [Double] -> [Double] -> [Double] -> Double
expectationValue dx bra ket quantity = (*dx') $ sum $ zipWith (*) bra $ zipWith (*) ket quantity
    where   dx' = abs dx


normalize :: Double -> [Double] -> [Double]
normalize dx wf = map (\w -> w/n) wf
    where   n = sqrt $ expectationValue dx wf wf (repeat 1) 

hypervirial :: Int -> Double -> (Double -> Double) -> (Double -> Double) -> Energy -> [Double] -> ([Double] -> Double) -> Double
hypervirial l mu force pot e xs exp = factor * (fExp + mu*4*l'*eExp)
    where   l' = fromIntegral l
            factor = 2*mu^(2*l+1)*((fromIntegral $ factorial l)/(2*l'+1))^2
            fExp = exp (map (\x -> (force (x/mu)) / (x^(2*l))) xs)
            eExp = exp (map (\x -> (e - pot (x/mu)) / x^(2*l+1)) xs)

factorial :: Int -> Int
factorial 0 = 1
factorial n = n * factorial (n-1)


solve :: Config -> (Double, Double)
solve config =  let xs = [xmax, xmax+dx..xmin]
                    cs = coefficients dx (potential l mu v xs)
                    chi' = ordinate cs
                    eev = head $ drop n $ scan chi' (emin, emax) de
                    e = energy chi' eev
                    (Just eng) = e
                    wf = wavefunction dx cs eng
                    exp' = expectationValue dx wf wf
                    origin = hypervirial l mu f v (eng*mu/2) xs exp'
                in  (eng*mu/2, origin)
    where   mu = mass config
            (Range xmax xmin dx) = xrange config
            (Range emin emax de) = erange config
            (QNumber n l) = qnumber config
            v = pot config
            f = force config
