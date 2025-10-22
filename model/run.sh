echo "Building..."
g++ -O3 -std=c++17 *.cpp -o netfil

if [ $? -eq 0 ]; then
  echo "Build successful. Running..."
  ./netfil out.csv 2 AS 1.0 1500
else
  echo "Build failed."
fi