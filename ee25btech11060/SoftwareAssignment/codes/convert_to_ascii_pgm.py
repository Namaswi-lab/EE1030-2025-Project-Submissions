from PIL import Image

# Input JPEG
input_path = "../figs/einstein.jpg"
# Output ASCII PGM
output_path = "../figs/einstein_ascii.pgm"

# Open image and convert to grayscale
img = Image.open(input_path).convert("L")
pixels = list(img.getdata())
width, height = img.size

# Write ASCII PGM manually
with open(output_path, "w") as f:
    f.write("P2\n")
    f.write(f"{width} {height}\n")
    f.write("255\n")
    for i, p in enumerate(pixels):
        f.write(f"{p} ")
        if (i + 1) % width == 0:
            f.write("\n")

print("✅ Successfully wrote ASCII PGM (P2) to:", output_path)

