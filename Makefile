
build:
	CC=/usr/bin/musl-gcc  cargo build --release --target=x86_64-unknown-linux-musl
#cross build --target x86_64-unknown-linux-musl

h:build
	target/x86_64-unknown-linux-musl/release/excord-lr --help

t:build 
	time target/x86_64-unknown-linux-musl/release/excord-lr -b test/mini_aln.sort.bam -o test/t 


tnomerge:build 
	time target/x86_64-unknown-linux-musl/release/excord-lr -n -b test/mini_aln.sort.bam -o test/t-nomerge 


te:
	time samtools view -b  test/mini_aln.sort.bam |  test/excord  --discordantdistance 500 /dev/stdin  >test/te

t2:build
	time target/x86_64-unknown-linux-musl/release/excord-lr -b test/full_aln.sort.bam -o test/tt 

tteclean:
	rm test/ttex

tte: 
	time samtools view -b  test/mini_aln.sort.bam |  \
	/mnt/e/projects/git_repo/excord/cmd/excord/excord \
	--discordantdistance 500 /dev/stdin  >test/ttex
ttes:
	cat test/ttex test/t | sort | grep -v 'align'  | grep -v 'align' | sort -k 1,1 -k 2,2 -k 5,5 -k 6,6 |tt |less -SN 


tte2: 
	time samtools view -b  test/mini_aln.sort.bam |  \
	test/excord \
	--discordantdistance 500 /dev/stdin  >test/ttex2
tte2s:
	cat test/ttex2 test/t | sort | grep -v 'align'  | grep -v 'align' | sort -k 1,1 -k 2,2 -k 5,5 -k 6,6 |tt |less -SN 

t2f:
	sudo /home/zhengxc/.cargo/bin/flamegraph  \
	-o my_flamegraph.svg  \
	-- \
	target/x86_64-unknown-linux-musl/release/excord-lr   \
	-b test/full_aln.sort.bam \
	-o test/tt 