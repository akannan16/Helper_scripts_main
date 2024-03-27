docker run \
 --rm \
 -it \
 -v /:/host \
 --mount type=bind,source=/folder1/,target=/folder2/ \
 path/to/docker/release/ver
bash
