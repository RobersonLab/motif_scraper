from setuptools import setup

setup(
	name = "motif_scraper",
	packages = ['motif_scraper'],
	version = "1.0.0",
	description = 'Tool for finding degenerate motifs in FASTA files',
	author = "Elisha Roberson",
	author_email = 'dr.eli.roberson@gmail.com',
	url = 'https://github.com/RobersonLab/motif_scraper',
	license = 'MIT',
	install_requires = ['pyfaidx', 'regex', 'six'],
	entry_points = {'console_scripts':["motif_scraper = motif_scraper.__main__:main"]},
	test_suite = 'nose.collector',
	tests_require = ['nose']
)
