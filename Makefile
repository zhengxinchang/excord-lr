
build:
	CC=/usr/bin/musl-gcc  cargo build --release --target=x86_64-unknown-linux-musl
#cross build --target x86_64-unknown-linux-musl

h:build
	target/x86_64-unknown-linux-musl/release/excord-lr --help

t:build 
	time target/x86_64-unknown-linux-musl/release/excord-lr -b test/mini_aln.sort.bam -o test/t 

te:
	time samtools view -b  test/mini_aln.sort.bam |  test/excord  --discordantdistance 500 /dev/stdin  >test/te

t2:build
	time target/x86_64-unknown-linux-musl/release/excord-lr -b test/full_aln.sort.bam -o test/tt 

tte:
	time samtools view -b  test/mini_aln.sort.bam |  test/excord  --discordantdistance 500 /dev/stdin  >test/ttex

t2f:
	sudo /home/zhengxc/.cargo/bin/flamegraph  \
	-o my_flamegraph.svg  \
	-- \
	target/x86_64-unknown-linux-musl/release/excord-lr   \
	-b test/full_aln.sort.bam \
	-o test/tt 