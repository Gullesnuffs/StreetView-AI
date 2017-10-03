from geopy.distance import vincenty

N, M, T, C, S, = map(int, input().split())

cords = [tuple(map(float,input().split())) for i in range(N)]
roads = [tuple(map(int,input().split())) for i in range(M)]

for c in cords:
  print(c)
for r in roads:
  A, B, D = r
  L = int(vincenty(cords[r[0]], cords[r[1]]).meters)
  # lets go with a speed of 15 km/h
  C = int(L / (1500 / 360))
  print(A, B, D, C, L)