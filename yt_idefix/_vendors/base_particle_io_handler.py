from typing import DefaultDict, List

from yt.utilities.io_handler import BaseIOHandler


# this class is vendored from yt 4.1dev
class BaseParticleIOHandler(BaseIOHandler):
    def _sorted_chunk_iterator(self, chunks):
        chunks = list(chunks)
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        yield from sorted(data_files, key=lambda x: (x.filename, x.start))

    def _count_particles_chunks(
        self,
        psize: DefaultDict[str, int],
        chunks,
        ptf: DefaultDict[str, List[str]],
        selector,
    ) -> DefaultDict[str, int]:
        if getattr(selector, "is_all_data", False):
            for data_file in self._sorted_chunk_iterator(chunks):
                for ptype in ptf.keys():
                    psize[ptype] += data_file.total_particles[ptype]
        else:
            # we must apply the selector and count the result
            psize = self._count_selected_particles(psize, chunks, ptf, selector)
        return psize
