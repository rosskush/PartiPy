import os

files = os.listdir()

# print(files)
gifs = [item for item in files if item.endswith('.gif')]

print(gifs)


for gif in gifs:
	mp4 = gif.replace('.gif','.mp4')
	if not os.path.exists(mp4):
		print(gif)
		arg = f'ffmpeg -i {gif} -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" {mp4}'
		os.system(arg)
		