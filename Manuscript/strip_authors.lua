-- strip_authors.lua: when rendered with -M anonymous=true, drop every author
-- and affiliation field so the document carries no author block.
function Meta(meta)
  if meta.anonymous == true or pandoc.utils.stringify(meta.anonymous or "") == "true" then
    for _, k in ipairs({"author", "authors", "by-author", "affiliations", "by-affiliation", "institute"}) do
      meta[k] = nil
    end
  end
  return meta
end
